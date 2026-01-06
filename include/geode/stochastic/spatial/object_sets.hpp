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

#pragma once

#include <absl/container/flat_hash_map.h>

#include <geode/basic/uuid.hpp>

#include <geode/geometry/basic_objects/segment.hpp>
#include <geode/geometry/bounding_box.hpp>
#include <geode/geometry/point.hpp>

#include <geode/stochastic/spatial/object_neighborhood.hpp>
#include <geode/stochastic/spatial/object_set.hpp>

namespace geode {
template <typename Type> struct ObjectRef {
  const Type &object;
  uuid set_id;
};
} // namespace geode

namespace geode {
template <typename Type> class ObjectSets {
  OPENGEODE_DISABLE_COPY(ObjectSets);

public:
  ObjectSets() noexcept = default;
  ObjectSets(ObjectSets &&) noexcept = default;
  ObjectSets &operator=(ObjectSets &&) noexcept = default;

  const ObjectSet<Type> &get_set(const uuid &set_id) const;
  const Type &get_object(const ObjectId &object_id) const;
  std::vector<ObjectId> get_all_object() const;
  std::vector<ObjectId> get_objects_in_set(const uuid &set_id) const;

  index_t nb_sets() const;
  index_t nb_objects_in_set(const uuid &set_id) const;
  index_t nb_objects() const;

  uuid add_set(std::string_view name);
  ObjectId add_object(Type &&object, const uuid &set_id, bool fixed);
  void update_free_object(const ObjectId &object_id, Type &&object);
  void remove_free_object(const ObjectId &object_id);

  // Object neighbor search by ObjectId (always excludes self)
  std::vector<ObjectId> neighbors(const ObjectId &object_id,
                                  const std::vector<uuid> &targeted_set_ids,
                                  double searching_distance) const;
  // Object neighbor search by arbitrary object (return self if in the
  // object_set)
  std::vector<ObjectId> neighbors(const Type &object,
                                  const std::vector<uuid> &targeted_set_ids,
                                  double searching_distance) const;

  std::string string() const;

private:
  ObjectSet<Type> &get_set(const uuid &set_id);
  //  void update_neighborhood_remove_context( const ObjectId& object_id
  //  );

private:
  absl::flat_hash_map<uuid, ObjectSet<Type>> sets_;
  ObjectNeighborhood<Type::dim> neighborhood_;
};
} // namespace geode
