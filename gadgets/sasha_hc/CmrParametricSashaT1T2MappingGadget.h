/**
\file   CmrParametricSashaT1T2MappingGadget.h
\brief  This is the class gadget for cardiac T1/T2 SASHA mapping, working on the IsmrmrdImageArray.
        If the T1p parameters are set, it will perform T1/T2/T1p mapping
\author Hui Xue
*/

#pragma once

#include "sashahclib_export.h"
#include "CmrParametricMappingGadget.h"

#ifndef GADGETRON_IMAGE_T1RHOMAP
    #define GADGETRON_IMAGE_T1RHOMAP                       "T1RHO"
#endif // GADGETRON_IMAGE_T1RHOMAP

namespace Gadgetron {

    class EXPORTSASHAHC CmrParametricSashaT1T2MappingGadget : public CmrParametricMappingGadget
    {
    public:
        GADGET_DECLARE(CmrParametricSashaT1T2MappingGadget);

        typedef CmrParametricMappingGadget BaseClass;

        typedef float real_value_type;
        typedef std::complex<real_value_type> ValueType;
        typedef Gadgetron::hoNDImage<real_value_type, 2> ImageType;
        typedef ValueType T;

        CmrParametricSashaT1T2MappingGadget();
        ~CmrParametricSashaT1T2MappingGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the mapping
        /// ------------------------------------------------------------------------------------

        GADGET_PROPERTY(max_iter,   size_t, "Maximal number of iterations",                   150);
        GADGET_PROPERTY(thres_func, double, "Threshold for minimal change of cost function", 1e-4);
        GADGET_PROPERTY(max_T1,     double, "Maximal T1 allowed in mapping (ms)",            4000);
        GADGET_PROPERTY(max_T2,     double, "Maximal T2 allowed in mapping (ms)",            1000);
        GADGET_PROPERTY(max_T1p,    double, "Maximal T1p allowed in mapping (ms)",           1000);

        GADGET_PROPERTY(anchor_image_index, size_t, "Index for anchor image; by default, the first image is the anchor (without SR pulse)", 0);
        GADGET_PROPERTY(anchor_TS,          double, "Saturation time for anchor",            100000);

        GADGET_PROPERTY(color_lut_t1map_15T, std::string, "Color lookup table for t1 map at 1.5T", "GadgetronT1_SR_1_5T.pal");
        GADGET_PROPERTY(color_lut_t1map_3T, std::string, "Color lookup table for t1 map at 1.5T", "GadgetronT1_SR_3T.pal");

        GADGET_PROPERTY(color_lut_t2map_15T, std::string, "Color lookup table for t2 map at 1.5T", "GadgetronT2_1_5T.pal");
        GADGET_PROPERTY(color_lut_t2map_3T, std::string, "Color lookup table for t2 map at 1.5T", "GadgetronT2_3T.pal");

        GADGET_PROPERTY(color_lut_t1pmap_15T, std::string, "Color lookup table for t1p map at 1.5T", "GadgetronT2_1_5T.pal");
        GADGET_PROPERTY(color_lut_t1pmap_3T, std::string, "Color lookup table for t1p map at 1.5T", "GadgetronT2_3T.pal");

        GADGET_PROPERTY(window_center_t1map_15T, double, "Window center for T1 map at 1.5T", 1300);
        GADGET_PROPERTY(window_width_t1map_15T,  double, "Window width for T1 map at 1.5T",  1300);
        GADGET_PROPERTY(window_center_t2map_15T, double, "Window center for T2 map at 1.5T",   60);
        GADGET_PROPERTY(window_width_t2map_15T,  double, "Window width for T2 map at 1.5T",   120);
        GADGET_PROPERTY(window_center_t1pmap_15T, double, "Window center for T1p map at 1.5T", 60);
        GADGET_PROPERTY(window_width_t1pmap_15T, double, "Window width for T1p map at 1.5T", 120);

        GADGET_PROPERTY(window_center_t1map_3T, double, "Window center for T1 map at 3T",       1250);
        GADGET_PROPERTY(window_width_t1map_3T,  double, "Window width for T1 map at 3T",        2500);
        GADGET_PROPERTY(window_center_t2map_3T, double, "Window center for T2 map at 3T",       60);
        GADGET_PROPERTY(window_width_t2map_3T,  double, "Window width for T2 map at 3T",        120);
        GADGET_PROPERTY(window_center_t1pmap_3T, double, "Window center for T1p map at 3T",     60);
        GADGET_PROPERTY(window_width_t1pmap_3T, double, "Window width for T1p map at 3T",       120);

        GADGET_PROPERTY(scaling_factor_t1map, double, "Scale factor for t1map", 1.0);
        GADGET_PROPERTY(scaling_factor_t2map, double, "Scale factor for t2map", 10.0);
        GADGET_PROPERTY(scaling_factor_t1pmap, double, "Scale factor for t1pmap", 10.0);

        GADGET_PROPERTY(has_HC, bool, "Whether to has HC lines", true);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // saturation recovery / T2 / T1rho prep times
        std::vector<float> prep_times_ts_;
        std::vector<float> prep_times_t2p_;
        std::vector<float> prep_times_t1p_;

        float              time_t2p_to_center_kspace_;
        std::vector<float> t2p_rf_duration_;

        size_t             num_rep_;

        bool               has_t1p_mapping_;

        // --------------------------------------------------
        // functional functions
        // --------------------------------------------------

        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1);

        virtual int close(unsigned long flags);

        // function to perform the mapping
        // data: input image array [RO E1 E2 CHA N S SLC]
        // map and map_sd: mapping result and its sd
        // para and para_sd: other parameters of mapping and its sd
        virtual int perform_mapping(IsmrmrdImageArray& data, IsmrmrdImageArray& map, IsmrmrdImageArray& para, IsmrmrdImageArray& map_sd, IsmrmrdImageArray& para_sd);

        virtual int perform_multi_mapping(IsmrmrdImageArray& data, IsmrmrdImageArray& t1map, IsmrmrdImageArray& t2map, IsmrmrdImageArray& para, IsmrmrdImageArray& map_sd, IsmrmrdImageArray& para_sd);
        virtual int perform_multi_mapping(IsmrmrdImageArray& data, IsmrmrdImageArray& t1map, IsmrmrdImageArray& t2map, IsmrmrdImageArray& t1pmap, IsmrmrdImageArray& para, IsmrmrdImageArray& map_sd, IsmrmrdImageArray& para_sd);
        virtual int fill_multi_map_header(IsmrmrdImageArray& t1map, IsmrmrdImageArray& t2map, IsmrmrdImageArray& t1pmap);

    };
}
