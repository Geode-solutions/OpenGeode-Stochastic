/*
 * Copyright (c) 2019 - 2025 Geode-solutions
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#include <geode/geometry/point.hpp>
#include <geode/stochastic/sampling/direct/object_set_sampler/point_set_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/helpers/simulation_runner.hpp>
#include <geode/stochastic/sampling/mcmc/metropolis_hasting_sampler.hpp>
#include <geode/stochastic/sampling/mcmc/models/components/density_term.hpp>
#include <geode/stochastic/sampling/mcmc/models/gibbs_energy.hpp>
#include <geode/stochastic/sampling/mcmc/proposal/classical_proposals.hpp>
#include <geode/stochastic/spatial/object_sets.hpp>
namespace {
struct SetDescription {
  std::string name;
  double birth_ratio{1.0};
  double death_ratio{1.0};
  double change_ratio{1.0};
};

struct PoissonDensityDescription {
  std::string name;
  double density;
  double target_count;
};

class PoissonSimulationRunner : public geode::SimulationRunner<geode::Point2D> {
public:
  PoissonSimulationRunner(const geode::SpatialDomain<2> &domain)
      : geode::SimulationRunner<geode::Point2D>(domain) {};

  void add_set_descriptor(const SetDescription &descriptor) {
    set_descriptors_.push_back(descriptor);
  }

  void add_density_descriptor(const PoissonDensityDescription &descriptor) {
    density_descriptors_.push_back(descriptor);
  }

  void initialize() override {
    auto proposal_kernel =
        std::make_unique<geode::ProposalKernel<geode::Point2D>>();

    // Mapping set names -> UUID
    std::unordered_map<std::string, geode::uuid> name_to_uuid;

    // Step 1: create object sets and samplers
    for (const auto &set_desc : set_descriptors_) {
      const auto set_id = this->object_sets_.add_set(set_desc.name);
      name_to_uuid[set_desc.name] = set_id;

      this->set_samplers_.push_back(
          std::make_unique<geode::UniformPointSetSampler<2>>(domain_));

      geode::add_birth_death_change_moves(
          this->set_samplers_.back(), *proposal_kernel, set_id,
          set_desc.birth_ratio, set_desc.death_ratio, set_desc.change_ratio);
    }

    // Step 2: create energy terms
    for (const auto &energy_desc : density_descriptors_) {
      const auto set_id = name_to_uuid.at(energy_desc.name);

      this->ordered_energy_terms_.push_back(
          this->energy_terms_collection_.add_energy_term(
              std::make_unique<geode::DensityTerm<geode::Point2D>>(
                  absl::StrCat(energy_desc.name, "_density"),
                  energy_desc.density, std::vector<geode::uuid>{set_id},
                  this->domain_)));

      this->ordered_target_statistics_.push_back(energy_desc.target_count);
    }

    this->mh_sampler_ =
        std::make_unique<geode::MetropolisHastings<geode::Point2D>>(
            this->energy_terms_collection_, std::move(proposal_kernel));
  }

  void
  check_statistics(const geode::StatisticsMonitor &statistic_monitoring) const {
    const auto &computed_means = statistic_monitoring.means();
    const auto &computed_variances = statistic_monitoring.variances();

    for (const auto stat_id :
         geode::Range{this->energy_terms_collection_.size()}) {
      const auto &term =
          energy_terms_collection_.get(ordered_energy_terms_[stat_id]);

      const auto expected_means = this->ordered_target_statistics_[stat_id];

      const auto target_vs_mean_error =
          std::abs(computed_means[stat_id] - expected_means) / expected_means;

      OPENGEODE_EXCEPTION(
          target_vs_mean_error < 0.05, "[MH test] statistic value ",
          computed_means[stat_id], " for energy term: ", term.name().data(),
          " not close enough to expected value ", expected_means,
          " --> error : ", target_vs_mean_error);

      const auto target_vs_variance_error =
          std::abs(computed_variances[stat_id] - expected_means) /
          expected_means;

      OPENGEODE_EXCEPTION(
          target_vs_variance_error < 0.15, "[MH test] variance of statistic ",
          computed_variances[stat_id], " for energy term: ", term.name().data(),
          " not close enough to expected value ", expected_means,
          " --> error : ", target_vs_variance_error);
    }
  }

private:
  geode::BoundingBox2D box_;
  std::vector<SetDescription> set_descriptors_;
  std::vector<PoissonDensityDescription> density_descriptors_;
};

void test_single_type_poisson() {
  geode::Logger::info("TEST - MH SINGLE TYPE POISSON");

  geode::RandomEngine engine;
  engine.set_seed("@mh-test-single-POISSON@");

  geode::BoundingBox2D box;
  box.add_point(geode::Point2D{{0.0, 0.0}});
  box.add_point(geode::Point2D{{10.0, 10.0}});
  geode::SpatialDomain domain(box, 0.);

  std::array<double, 4> birth_ratio{0.1, 0.5, 2., 4.};
  std::array<double, 4> change_ratio{0., 1., 1., 0.};

  for (const auto config : geode::Range{birth_ratio.size()}) {
    // --- Set description
    SetDescription setA;
    setA.name = "A";
    setA.birth_ratio = birth_ratio[config];
    setA.death_ratio = 1.0;
    setA.change_ratio = change_ratio[config];

    // --- Energy term description
    PoissonDensityDescription densityA;
    densityA.name = "A";
    densityA.density = 0.3;
    densityA.target_count = 30.0;

    PoissonSimulationRunner runner(domain);
    runner.add_set_descriptor(setA);
    runner.add_density_descriptor(densityA);
    runner.initialize();

    // run simulation
    geode::SimulationPrinterConfigurator printer_config;
    printer_config.output_folder = absl::StrCat(
        printer_config.output_folder, "/sim_point_poisson_test_", config);

    geode::SimulationConfigurator sim_config;
    sim_config.realizations = 1000;
    sim_config.metropolis_hasting_steps = 1000;
    sim_config.burn_in_steps = 1000;
    sim_config.printer = printer_config;

    auto statistic_monitoring = runner.run(engine, sim_config);
    runner.check_statistics(statistic_monitoring);
  }

  geode::Logger::info("--> SUCCESS!");
}

void test_multitype_poisson() {
  geode::Logger::info("TEST - MH MULTITYPE POISSON");

  geode::RandomEngine engine;
  engine.set_seed("@mh-test-POISSON-multi@");

  geode::BoundingBox2D box;
  box.add_point(geode::Point2D{{0.0, 0.0}});
  box.add_point(geode::Point2D{{10.0, 10.0}});
  geode::SpatialDomain domain(box, 0.);

  // --- Set descriptions
  SetDescription set01{"set01", 2.0, 3.0, 1.0};
  SetDescription set02{"set02", 3.0, 0.5, 1.0};
  SetDescription set03{"set03", 4.0, 1.0, 1.0};

  // --- Energy term descriptions
  PoissonDensityDescription density01{"set01", 0.1, 10.0};
  PoissonDensityDescription density02{"set02", 0.4, 40.0};
  PoissonDensityDescription density03{"set03", 0.3, 30.0};

  PoissonSimulationRunner runner(domain);
  runner.add_set_descriptor(set01);
  runner.add_set_descriptor(set02);
  runner.add_set_descriptor(set03);

  runner.add_density_descriptor(density01);
  runner.add_density_descriptor(density02);
  runner.add_density_descriptor(density03);

  runner.initialize();

  // run simulation
  geode::SimulationPrinterConfigurator printer_config;
  printer_config.output_folder = absl::StrCat(
      printer_config.output_folder, "/sim_point_multitype_poisson_test");

  geode::SimulationConfigurator sim_config;
  sim_config.realizations = 1000;
  sim_config.metropolis_hasting_steps = 1000;
  sim_config.burn_in_steps = 1000;
  sim_config.printer = printer_config;

  auto statistic_monitoring = runner.run(engine, sim_config);
  runner.check_statistics(statistic_monitoring);

  geode::Logger::info("--> SUCCESS!");
}
} // namespace

int main() {
  try {
    geode::StochasticLibrary::initialize();
    geode::Logger::set_level(geode::Logger::LEVEL::debug);
    test_single_type_poisson();
    test_multitype_poisson();
    return 0;
  } catch (...) {
    return geode::geode_lippincott();
  }
}