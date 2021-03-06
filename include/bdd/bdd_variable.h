#pragma once
#include <array>

namespace LPMP {

    template<typename DERIVED>
    class bdd_variable_base {
        public:
            std::size_t first_node_index = std::numeric_limits<std::size_t>::max();
            std::size_t last_node_index = std::numeric_limits<std::size_t>::max();

            DERIVED* prev = nullptr;
            DERIVED* next = nullptr;

            std::size_t nr_bdd_nodes() const { return last_node_index - first_node_index; }
            bool is_first_bdd_variable() const { return prev == nullptr; }
            bool is_last_bdd_variable() const { return next == nullptr; }
            // bool is_initial_state() const { return *this == bdd_variable<DERIVED>{}; }
            friend bool operator==(const bdd_variable_base<DERIVED>&, const bdd_variable_base<DERIVED>&);
    };

    template<typename DERIVED>
    bool operator==(const bdd_variable_base<DERIVED>& x, const bdd_variable_base<DERIVED>& y)
    {
        return (x.first_node_index == y.first_node_index &&
            x.last_node_index == y.last_node_index &&
            x.prev == y.prev &&
            x.next == y.next); 
    }

    class bdd_variable : public bdd_variable_base<bdd_variable> {};

    template<typename DERIVED>
    class bdd_variable_mma_base : public bdd_variable_base<DERIVED> {
        public:
            double cost = std::numeric_limits<double>::infinity();
    };

    template<typename DERIVED>
    bool operator==(const bdd_variable_mma_base<DERIVED>& x, const bdd_variable_mma_base<DERIVED>& y)
    {
        return (x.first_node_index == y.first_node_index &&
            x.last_node_index == y.last_node_index &&
            x.prev == y.prev &&
            x.next == y.next &&
            x.cost == y.cost); 
    }

    class bdd_variable_mma : public bdd_variable_mma_base<bdd_variable_mma> {   
    };

    class bdd_variable_fix : public bdd_variable_mma_base<bdd_variable_fix> {
        public:
            std::size_t nr_feasible_low_arcs;
            std::size_t nr_feasible_high_arcs;
            std::size_t variable_index;
    };

    bool operator==(const bdd_variable_fix& x, const bdd_variable_fix& y)
    {
        return (x.first_node_index == y.first_node_index &&
            x.last_node_index == y.last_node_index &&
            x.prev == y.prev &&
            x.next == y.next &&
            x.cost == y.cost &&
            x.nr_feasible_low_arcs == y.nr_feasible_low_arcs &&
            x.nr_feasible_high_arcs == y.nr_feasible_high_arcs &&
            x.variable_index == y.variable_index); 
    }
}

