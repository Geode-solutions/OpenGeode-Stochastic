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

#include <geode/basic/identifier.hpp>
#include <geode/stochastic/common.hpp>

namespace geode
{
    template < typename Type >
    class ObjectSet : public Identifier
    {
    public:
        ObjectSet() noexcept = default;
        ObjectSet( ObjectSet&& ) noexcept = default;
        ObjectSet& operator=( ObjectSet&& ) noexcept = default;

        void set_name( std::string_view name );
        const Type& get_object( index_t index ) const;

        index_t add_fixed_object( Type&& object );
        void remove_fixed_object( index_t index );

        index_t add_free_object( Type&& object );
        void update_free_object( index_t index, Type&& object );
        void remove_free_object( index_t index );

        index_t nb_objects() const;
        index_t nb_fixed_objects() const;
        index_t nb_free_objects() const;

        bool empty() const;
        bool is_fixed( index_t index ) const;

        std::string string() const;

    private:
        void do_remove_object( index_t index );

    private:
        std::vector< Type > objects_;
        index_t first_free_object_{ 0 };
    };
} // namespace geode
