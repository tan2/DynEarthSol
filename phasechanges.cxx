#include "constants.hpp"
#include "parameters.hpp"
#include "markerset.hpp"
#include "utils.hpp"

#include "phasechanges.hpp"

namespace {

    // A template to implement phase_change_fn
    int custom_phase_changes(const Param& param, const Variables& var,
                             const MarkerSet& ms, int m)
    {
        int current_mt = ms.get_mattype(m);
        int e = ms.get_elem(m);
        const double* eta = ms.get_eta(m);

        // Get temperature at the marker
        double T = 0;
        const int *conn = (*var.connectivity)[e];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            T += (*var.temperature)[conn[i]] * eta[i];
        }

        // Get pressure, which is constant in the element
        double P = trace((*var.stress)[e]) / NDIMS;

        // Set new mattype the same as current mattype for now
        int new_mt = current_mt;

        switch (current_mt) {
        case 0:
            break;
        case 1:
            break;
        }

        return new_mt;
    }

} // anonymous namespace


void phase_changes(const Param& param, const Variables& var,
                   MarkerSet& ms, int_vec2D& elemmarkers)
{
    if (param.mat.nmat == 1 || param.mat.phase_change_option == 0) return;

    int (*phase_change_fn)(const Param&, const Variables&, const MarkerSet&, int) = NULL;
    switch (param.mat.phase_change_option) {
    case 101:
        phase_change_fn = custom_phase_changes;
        break;
    default:
        std::cerr << "Error: unknown phase_change_option: " << param.mat.phase_change_option << '\n';
        std::exit(1);
    }

    #pragma omp parallel for default(none)          \
        shared(param, var, ms, elemmarkers, phase_change_fn)
    for (int m=0; m<ms.get_nmarkers(); m++) {
        int current_mt = ms.get_mattype(m);
        int new_mt = phase_change_fn(param, var, ms, m);

        if (new_mt != current_mt) {
            int e = ms.get_elem(m);
            ms.set_mattype(m, new_mt);
            // prevent concurrent modification on elemmarkers
            #pragma omp critical (phase_changes)
            {
                --elemmarkers[e][current_mt];
                ++elemmarkers[e][new_mt];
            }
        }

    }
}
