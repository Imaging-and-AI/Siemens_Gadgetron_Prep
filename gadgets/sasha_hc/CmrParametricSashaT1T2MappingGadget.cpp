
#include "CmrParametricSashaT1T2MappingGadget.h"
#include <iomanip>
#include <sstream>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_utility.h"
// #include "cmr_t1_mapping.h"
#include "cmr_sasha_t1_t2_mapping.h"

namespace Gadgetron {

    CmrParametricSashaT1T2MappingGadget::CmrParametricSashaT1T2MappingGadget() : BaseClass()
    {
        this->has_t1p_mapping_ = false;
    }

    CmrParametricSashaT1T2MappingGadget::~CmrParametricSashaT1T2MappingGadget()
    {
    }

    int CmrParametricSashaT1T2MappingGadget::process_config(ACE_Message_Block* mb)
    {
        GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget::process_config(...) starts ... ");
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        GDEBUG_STREAM("meas_max_idx_.repetition is " << meas_max_idx_.repetition);
        GDEBUG_STREAM("meas_max_idx_.average is " << meas_max_idx_.average);
        GDEBUG_STREAM("meas_max_idx_.set is " << meas_max_idx_.set);
        num_encoding_spaces_ = 1;//h.encoding.size();

        num_rep_ = meas_max_idx_.repetition + 1;
        if(num_rep_==1 && meas_max_idx_.average>=1)
        {
            num_rep_ = meas_max_idx_.average + 1;
        }
        GDEBUG_STREAM("num_rep_ is " << num_rep_);

        // Each image in the set must have a prep time
        // this->prep_times_ts_.resize( this->meas_max_idx_.set + 1);
        // this->prep_times_t2p_.resize(this->meas_max_idx_.set + 1);
        this->prep_times_ts_.resize(  (this->meas_max_idx_.set + 1) * (num_rep_));
        this->prep_times_t2p_.resize( (this->meas_max_idx_.set + 1) * (num_rep_));
        this->prep_times_t1p_.resize( (this->meas_max_idx_.set + 1) * (num_rep_));
        this->t2p_rf_duration_.resize((this->meas_max_idx_.set + 1) * (num_rep_));

        this->time_t2p_to_center_kspace_ = 0;
        float t2p_rf_duration_old = 0.0; // Backwards compatibility with old data with only a single T2p duration

        //if (this->imaging_prep_time_from_protocol.value())
        {
            // Read in minimum time to center (from start of single-shot readout to center k-space)
            if (h.userParameters->userParameterDouble.size() > 0)
            {
                std::vector<ISMRMRD::UserParameterDouble>::const_iterator iter = h.userParameters->userParameterDouble.begin();

                for (; iter != h.userParameters->userParameterDouble.end(); iter++)
                {
                    std::string usrParaName  = iter->name;
                    long        usrParaValue = iter->value;

                    std::string strTimeToCenter("TimeT2pToCenterKspace");
                    if (usrParaName == strTimeToCenter)
                    {
                        this->time_t2p_to_center_kspace_ = (float)usrParaValue;
                        GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, found time T2p to center : " << this->time_t2p_to_center_kspace_);
                    }

                    std::string strT2pRfDuration("T2pRfDuration");
                    if (usrParaName == strT2pRfDuration)
                    {
                        t2p_rf_duration_old = (float)usrParaValue;
                        GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, found single T2p RF duration (old format) : " << t2p_rf_duration_old);
                    }
                }
            }

            // First SASHA image is an anchor unless otherwise specified
            this->prep_times_ts_[ 0] = 10000000; // 10 s dummy value
            this->prep_times_t2p_[0] = 0;
            this->prep_times_t1p_[0] = 0;

            // Read in T1 and T2 prep times
            // FIXME: This doesn't allow for the first image to be T2p or SR
            size_t iT1    = 1;
            size_t iT2    = 1;
            size_t iT1p   = 1;
            size_t iT2pRf = 1;
            if (h.userParameters->userParameterDouble.size() > 0)
            {
                std::vector<ISMRMRD::UserParameterDouble>::const_iterator iter = h.userParameters->userParameterDouble.begin();

                for (; iter != h.userParameters->userParameterDouble.end(); iter++)
                {
                    std::string usrParaName  = iter->name;
                    double      usrParaValue = iter->value;

                    std::stringstream strT1, strT2, strT1p, strT2pRfDuration;
                    strT1            << "SatRecTime_"      << iT1;
                    strT2            << "T2PrepDuration_"  << iT2;
                    strT1p           << "T1pPrepDuration_" << iT1p;
                    strT2pRfDuration << "T2pRfDuration_"   << iT2pRf;

                    // GDEBUG_STREAM("Searching for " << strT1.str() << " and "<< strT2.str() << " and " << strT1p.str());

                    if (usrParaName == strT1.str() && iT1 <= this->meas_max_idx_.set)
                    {
                        GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, found SR prep time : " << iT1 << " - " << usrParaValue);
                        for (size_t i = 0; i < num_rep_; i++)
                        {
                            size_t ind = iT1+i*(this->meas_max_idx_.set+1);
                            if (ind < this->prep_times_ts_.size())
                            {
                                this->prep_times_ts_[iT1+i*(this->meas_max_idx_.set+1)] = (float)usrParaValue;
                            }
                        }
                        iT1++;
                    }
                    else if (usrParaName == strT2.str() && iT2 <= this->meas_max_idx_.set)
                    {
                        GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, found T2 prep time : " << iT2 << " - " << usrParaValue);
                        for (size_t i = 0; i < num_rep_; i++)
                        {
                            size_t ind = iT2+i*(this->meas_max_idx_.set+1);
                            if (ind < this->prep_times_t2p_.size())
                            {
                                this->prep_times_t2p_[iT2+i*(this->meas_max_idx_.set+1)] = (float)usrParaValue;
                            }
                        }
                        iT2++;
                    }
                    else if (usrParaName == strT1p.str() && iT1p <= this->meas_max_idx_.set)
                    {
                        GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, found T1p prep time : " << iT1p << " - " << usrParaValue);
                        for (size_t i = 0; i < num_rep_; i++)
                        {
                            size_t ind = iT1p + i * (this->meas_max_idx_.set + 1);
                            if (ind < this->prep_times_t1p_.size())
                            {
                                this->prep_times_t1p_[iT1p + i * (this->meas_max_idx_.set + 1)] = (float)usrParaValue;
                            }
                        }
                        iT1p++;
                    }
                    else if (usrParaName == strT2pRfDuration.str() && iT2pRf <= this->meas_max_idx_.set)
                    {
                        GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, found T2p RF duration : " << iT2pRf << " - " << usrParaValue);
                        for (size_t i = 0; i < num_rep_; i++)
                        {
                            size_t ind = iT2pRf + i * (this->meas_max_idx_.set + 1);
                            if (ind < this->t2p_rf_duration_.size())
                            {
                                this->t2p_rf_duration_[iT2pRf + i * (this->meas_max_idx_.set + 1)] = (float)usrParaValue;
                            }
                        }
                        iT2pRf++;
                    }
                }
            }

            // Populate t2p_rf_duration_ with old format data where only a single value was exported
            for (size_t i = 0; i < this->t2p_rf_duration_.size(); i++)
            {
                if (this->prep_times_t2p_[i] > 0)
                {
                    this->t2p_rf_duration_[i] = t2p_rf_duration_old;
                }
            }

            // First SASHA image is an anchor unless otherwise specified
            if (this->has_HC.value())
            {
                this->prep_times_ts_[ 0] = 10000000;
                this->prep_times_t2p_[0] = 0;
            }
            else
            {
                this->prep_times_ts_[ 0] = 10000000;

                this->prep_times_t2p_[0] = this->prep_times_t2p_[1];
                size_t SET = this->meas_max_idx_.set+1;
                for (size_t n=0; n<num_rep_; n++)
                {
                    this->prep_times_t2p_[n*SET] = this->prep_times_t2p_[n*SET+1];
                }
            }

            // Check if there are any images with T1rho prep
            for (size_t i = 0; i < this->prep_times_t1p_.size(); i++)
            {
                if (this->prep_times_t1p_[i] > 0)
                {
                    has_t1p_mapping_ = true;
                }
            }

            if (!has_t1p_mapping_)
            {
                GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, CANNOT find T1rho prep duration ... ");
            }
        }

        return GADGET_OK;
    }

    // Overloaded from CmrParametricMappingGadget in order to send both T1 and T2 data
    int CmrParametricSashaT1T2MappingGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("CmrParametricSashaT1T2MappingGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "CmrParametricSashaT1T2MappingGadget::process(...) starts ... ");

        // -------------------------------------------------------------

        process_called_times_++;

        // -------------------------------------------------------------

        IsmrmrdImageArray* data = m1->getObjectPtr();

        // print out data info
        if (verbose.value())
        {
            GDEBUG_STREAM("----> CmrParametricSashaT1T2MappingGadget::process(...) has been called " << process_called_times_ << " times ...");
            std::stringstream os;
            data->data_.print(os);
            GDEBUG_STREAM(os.str());
        }

        // -------------------------------------------------------------

        // some images do not need mapping
        if (data->meta_[0].length(skip_processing_meta_field.value().c_str())>0)
        {
            for (size_t ii = 0; ii < data->headers_.get_number_of_elements(); ii++)
            {
                data->headers_(ii).contrast = 0;
                GDEBUG_STREAM("ori image " << ii << ", pmu time is " << data->headers_(ii).physiology_time_stamp[0]);
            }

            if (this->next()->putq(m1) == -1)
            {
                GERROR("CmrParametricSashaT1T2MappingGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }
            GDEBUG_STREAM("Skipped image processing due to flag");
            return GADGET_OK;
        }

        // -------------------------------------------------------------

        size_t encoding = (size_t)data->meta_[0].as_long("encoding", 0);
        GDEBUG_STREAM("encoding:" << encoding);
        GDEBUG_STREAM("num_encoding_spaces_:" << num_encoding_spaces_);
        // GADGET_CHECK_RETURN(encoding<num_encoding_spaces_, GADGET_FAIL);

        std::string dataRole = std::string(data->meta_[0].as_str(GADGETRON_DATA_ROLE));

        size_t RO = data->data_.get_size(0);
        size_t E1 = data->data_.get_size(1);
        size_t E2 = data->data_.get_size(2);
        size_t CHA = data->data_.get_size(3);
        size_t N = data->data_.get_size(4);
        size_t S = data->data_.get_size(5);
        size_t SLC = data->data_.get_size(6);

        GDEBUG_STREAM("[RO E1 E2 CHA N S SLC] = " << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC);
        GDEBUG_STREAM("dataRole:" << dataRole);

        size_t ro, e1, e2, cha, n, s, slc;

        if (!this->has_HC.value())
        {
            hoNDArray<T> data_n_s;
            data_n_s.create(RO, E1, E2, CHA, N*S, 1, SLC);
            for (slc=0; slc<SLC; slc++)
            {
                for (s=0; s<S; s++)
                {
                    for (n=0; n<N; n++)
                    {
                        memcpy(&data_n_s(0, 0, 0, 0, n+s*N, 0, slc), &data->data_(0, 0, 0, 0, n, s, slc), sizeof(T)*RO*E1);
                    }
                }
            }

            data->data_ = data_n_s;
            N *= S;
            S = 1;
        }

        std::stringstream os;
        os << "_encoding_" << encoding << "_processing_" << process_called_times_;
        std::string str = os.str();

        // -------------------------------------------------------------

        if (!debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(data->data_, debug_folder_full_path_ + "data" + str);
        }

        // -------------------------------------------------------------
        // // if prep times are not read from protocol, find them from image header
        // if (!imaging_prep_time_from_protocol.value())
        // {
        //     this->prep_times_.resize(N, 0);

        //     size_t n;
        //     for (n = 0; n < N; n++)
        //     {
        //         this->prep_times_[n] = data->headers_(n).user_int[7] * 1e-3; // convert microsecond to ms
        //     }

        //     if (this->verbose.value())
        //     {
        //         GDEBUG_STREAM("Prep times are read from images ... ");
        //         for (n = 0; n < N; n++)
        //         {
        //             GDEBUG_STREAM("Image " << n << " - " << this->prep_times_[n] << " ms ");
        //         }
        //     }
        // }

        // Find overridden TS times from image header
        for (size_t n = 0; n < N; n++)
        {
            if (data->headers_(n).user_int[7] != 0)
            {
                this->prep_times_ts_[n] = data->headers_(n).user_int[7] * 1e-3; // convert microsecond to ms
                GDEBUG_STREAM("set ts from user_int, image "   <<                                       std::setw(2) << n
                           << ": TS = "  << std::fixed << std::setprecision(1) << std::setw(6) << this->prep_times_ts_[n]  << "* ms "
                           << ", T2p = " << std::fixed << std::setprecision(1) << std::setw(5) << this->prep_times_t2p_[n] << " ms "
                           << ", TSL = " << std::fixed << std::setprecision(1) << std::setw(5) << this->prep_times_t1p_[n] << " ms ");
            }
            else
            {
                GDEBUG_STREAM("Image "   <<                                       std::setw(2) << n
                           << ": TS = "  << std::fixed << std::setprecision(1) << std::setw(6) << this->prep_times_ts_[n]  << "  ms "
                           << ", T2p = " << std::fixed << std::setprecision(1) << std::setw(5) << this->prep_times_t2p_[n] << " ms "
                           << ", TSL = " << std::fixed << std::setprecision(1) << std::setw(5) << this->prep_times_t1p_[n] << " ms ");
            }
        }

        // -------------------------------------------------------------

        // calling the mapping
        // IsmrmrdImageArray map, para, map_sd, para_sd;
        IsmrmrdImageArray t1map, t2map, t1pmap, para, map_sd, para_sd;

        // int status = this->perform_mapping(*data, map, para, map_sd, para_sd);
        int status = 0;
        if(this->has_t1p_mapping_)
            status = this->perform_multi_mapping(*data, t1map, t2map, t1pmap, para, map_sd, para_sd);
        else
            status = this->perform_multi_mapping(*data, t1map, t2map, para, map_sd, para_sd);

        if (status != GADGET_OK)
        {
            GWARN("CmrParametricSashaT1T2MappingGadget::process, process incoming data failed ... ");

            for (size_t ii = 0; ii < data->headers_.get_number_of_elements(); ii++)
            {
                data->headers_(ii).contrast = 0;

                GDEBUG_STREAM("moco image " << ii << ", pmu time is " << data->headers_(ii).physiology_time_stamp[0]);
            }

            // sending the incoming images
            if (this->next()->putq(m1) == -1)
            {
                GERROR("CmrParametricSashaT1T2MappingGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }
        }
        else
        {
            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(t1map.data_,   debug_folder_full_path_ + "t1map"   + str);
                gt_exporter_.export_array_complex(t2map.data_,   debug_folder_full_path_ + "t2map"   + str);
                gt_exporter_.export_array_complex(para.data_,    debug_folder_full_path_ + "para"    + str);
                gt_exporter_.export_array_complex(map_sd.data_,  debug_folder_full_path_ + "map_sd"  + str);
                gt_exporter_.export_array_complex(para_sd.data_, debug_folder_full_path_ + "para_sd" + str);

                if (this->has_t1p_mapping_)
                {
                    gt_exporter_.export_array_complex(t1pmap.data_, debug_folder_full_path_ + "t1pmap" + str);
                }
            }

            // before sending the images, make sure the TS and TE are correct
            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        size_t ind = n + s * N + slc * N * S;

                        float ts = this->prep_times_ts_[n + s * N];
                        float t2p = this->prep_times_t2p_[n + s * N];

                        m1->getObjectPtr()->headers_(ind).contrast = 0;

                        m1->getObjectPtr()->meta_[ind].set(GADGETRON_IMAGE_SATURATIONTIME, (double)ts);
                        m1->getObjectPtr()->meta_[ind].set(GADGETRON_IMAGE_INVERSIONTIME, (double)ts);
                        m1->getObjectPtr()->meta_[ind].set(GADGETRON_IMAGE_ECHOTIME, (double)t2p);

                        if (this->has_t1p_mapping_)
                        {
                            float t1p = this->prep_times_t1p_[n + s * N];
                            m1->getObjectPtr()->meta_[ind].set("GADGETRON_T1RHO_PREP_TIME", (double)t1p);
                        }

                        GDEBUG_STREAM("moco image " << ind << ", pmu time is " << m1->getObjectPtr()->headers_(ind).physiology_time_stamp[0]);
                    }
                }
            }

            // sending the incoming images
            if (this->next()->putq(m1) == -1)
            {
                GERROR("CmrParametricSashaT1T2MappingGadget::process, passing incoming image array on to next gadget");
                return GADGET_FAIL;
            }

            GDEBUG_CONDITION_STREAM(verbose.value(), "CmrParametricSashaT1T2MappingGadget::process, deciding on whether to send maps");
            if (send_map.value() && (this->fill_multi_map_header(t1map, t2map, t1pmap) == GADGET_OK))
            {
                GDEBUG_CONDITION_STREAM(verbose.value(), "CmrParametricSashaT1T2MappingGadget::process, sending maps...");
                Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm1 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();
                *(cm1->getObjectPtr()) = t1map;
                if (this->next()->putq(cm1) == -1)
                {
                    GERROR("CmrParametricSashaT1T2MappingGadget::process, passing map on to next gadget");
                    return GADGET_FAIL;
                }

                Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm2 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();
                *(cm2->getObjectPtr()) = t2map;
                if (this->next()->putq(cm2) == -1)
                {
                    GERROR("CmrParametricSashaT1T2MappingGadget::process, passing map on to next gadget");
                    return GADGET_FAIL;
                }

                if (this->has_t1p_mapping_)
                {
                    Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm3 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();
                    *(cm3->getObjectPtr()) = t1pmap;
                    if (this->next()->putq(cm3) == -1)
                    {
                        GERROR("CmrParametricSashaT1T2MappingGadget::process, passing map on to next gadget");
                        return GADGET_FAIL;
                    }
                }
            }
        }

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    int CmrParametricSashaT1T2MappingGadget::fill_multi_map_header(IsmrmrdImageArray& t1map, IsmrmrdImageArray& t2map, IsmrmrdImageArray& t1pmap)
    {
        try
        {
            // Assume that the T1 and T2 maps have the same size
            size_t RO  = t1map.data_.get_size(0);
            size_t E1  = t1map.data_.get_size(1);
            size_t E2  = t1map.data_.get_size(2);
            size_t CHA = t1map.data_.get_size(3);
            size_t N   = t1map.data_.get_size(4);
            size_t S   = t1map.data_.get_size(5);
            size_t SLC = t1map.data_.get_size(6);

            size_t e2, cha, n, s, slc;

            std::string lut_t1map = color_lut_t1map_15T.value();
            std::string lut_t2map = color_lut_t2map_15T.value();
            std::string lut_t1pmap = color_lut_t1pmap_15T.value();
            double window_center_t1map = window_center_t1map_15T.value();
            double window_width_t1map  = window_width_t1map_15T.value();
            double window_center_t2map = window_center_t2map_15T.value();
            double window_width_t2map  = window_width_t2map_15T.value();
            double window_center_t1pmap = window_center_t1pmap_15T.value();
            double window_width_t1pmap = window_width_t1pmap_15T.value();

            if (this->field_strength_T_ > 2)
            {
                lut_t1map = color_lut_t1map_3T.value();
                lut_t2map = color_lut_t2map_3T.value();
                lut_t1pmap = color_lut_t1pmap_3T.value();
                window_center_t1map = window_center_t1map_3T.value();
                window_width_t1map  = window_width_t1map_3T.value();
                window_center_t2map = window_center_t2map_3T.value();
                window_width_t2map  = window_width_t2map_3T.value();
                window_center_t1pmap = window_center_t1pmap_3T.value();
                window_width_t1pmap = window_width_t1pmap_3T.value();
            }

            std::ostringstream ostr;
            ostr << "x" << (double)scaling_factor_t1map.value();
            std::string scalingStr_t1map = ostr.str();

            std::ostringstream ostr_unit;
            ostr_unit << std::setprecision(3) << 1.0f / scaling_factor_t1map.value() << "ms";
            std::string unitStr_t1map = ostr_unit.str();

            std::ostringstream ostr_t2map;
            ostr_t2map << "x" << (double)scaling_factor_t2map.value();
            std::string scalingStr_t2map = ostr_t2map.str();

            std::ostringstream ostr_unit_t2map;
            ostr_unit_t2map << std::setprecision(3) << 1.0f / scaling_factor_t2map.value() << "ms";
            std::string unitStr_t2map = ostr_unit_t2map.str();

            std::ostringstream ostr_t1pmap;
            ostr_t1pmap << "x" << (double)scaling_factor_t1pmap.value();
            std::string scalingStr_t1pmap = ostr_t1pmap.str();

            std::ostringstream ostr_unit_t1pmap;
            ostr_unit_t1pmap << std::setprecision(3) << 1.0f / scaling_factor_t1pmap.value() << "ms";
            std::string unitStr_t1pmap = ostr_unit_t1pmap.str();

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        size_t offset = n + s*N + slc*N*S;

                        t1map.headers_[offset].repetition = 0;
                        t1map.headers_[offset].average = 0;
                        t1map.headers_[offset].contrast = 0;

                        t1map.meta_[offset].set(GADGETRON_IMAGE_SCALE_RATIO,  (double)                     scaling_factor_t1map.value());
                        t1map.meta_[offset].set(GADGETRON_IMAGE_WINDOWCENTER, (long)(window_center_t1map * scaling_factor_t1map.value()));
                        t1map.meta_[offset].set(GADGETRON_IMAGE_WINDOWWIDTH,  (long)(window_width_t1map  * scaling_factor_t1map.value()));
                        t1map.meta_[offset].set(GADGETRON_IMAGE_COLORMAP,     lut_t1map.c_str());

                        t1map.meta_[offset].set(   GADGETRON_IMAGECOMMENT, t1map.meta_[offset].as_str(GADGETRON_DATA_ROLE));
                        t1map.meta_[offset].append(GADGETRON_IMAGECOMMENT, scalingStr_t1map.c_str());
                        t1map.meta_[offset].append(GADGETRON_IMAGECOMMENT, unitStr_t1map.c_str());

                        // t1map.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1MAP);
                        // t1map.meta_[offset].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1MAP);

                        GDEBUG_STREAM("T1 map, pmu time is " << t1map.headers_[offset].physiology_time_stamp[0]);

                        t2map.headers_[offset].repetition = 0;
                        t2map.headers_[offset].average = 0;
                        t2map.headers_[offset].contrast = 0;

                        t2map.meta_[offset].set(GADGETRON_IMAGE_SCALE_RATIO,  (double)                     scaling_factor_t2map.value());
                        t2map.meta_[offset].set(GADGETRON_IMAGE_WINDOWCENTER, (long)(window_center_t2map * scaling_factor_t2map.value()));
                        t2map.meta_[offset].set(GADGETRON_IMAGE_WINDOWWIDTH,  (long)(window_width_t2map  * scaling_factor_t2map.value()));
                        t2map.meta_[offset].set(GADGETRON_IMAGE_COLORMAP,     lut_t2map.c_str());

                        t2map.meta_[offset].set(   GADGETRON_IMAGECOMMENT, t2map.meta_[offset].as_str(GADGETRON_DATA_ROLE));
                        t2map.meta_[offset].append(GADGETRON_IMAGECOMMENT, scalingStr_t2map.c_str());
                        t2map.meta_[offset].append(GADGETRON_IMAGECOMMENT, unitStr_t2map.c_str());

                        // t2map.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T2MAP);
                        // t2map.meta_[offset].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T2MAP);

                        GDEBUG_STREAM("T2 map, pmu time is " << t2map.headers_[offset].physiology_time_stamp[0]);

                        if(this->has_t1p_mapping_)
                        {
                            t1pmap.headers_[offset].repetition = 0;
                            t1pmap.headers_[offset].average = 0;
                            t1pmap.headers_[offset].contrast = 0;

                            t1pmap.meta_[offset].set(GADGETRON_IMAGE_SCALE_RATIO, (double)scaling_factor_t1pmap.value());
                            t1pmap.meta_[offset].set(GADGETRON_IMAGE_WINDOWCENTER, (long)(window_center_t1pmap * scaling_factor_t1pmap.value()));
                            t1pmap.meta_[offset].set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(window_width_t1pmap * scaling_factor_t1pmap.value()));
                            t1pmap.meta_[offset].set(GADGETRON_IMAGE_COLORMAP, lut_t1pmap.c_str());

                            t1pmap.meta_[offset].set(GADGETRON_IMAGECOMMENT, t1pmap.meta_[offset].as_str(GADGETRON_DATA_ROLE));
                            t1pmap.meta_[offset].append(GADGETRON_IMAGECOMMENT, scalingStr_t1pmap.c_str());
                            t1pmap.meta_[offset].append(GADGETRON_IMAGECOMMENT, unitStr_t1pmap.c_str());

                            // we will need a new data role type
                            // t1pmap.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1RHOMAP);
                            // t1pmap.meta_[offset].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1RHOMAP);

                            t1pmap.meta_[offset].append(GADGETRON_IMAGE_INTENSITY_UNCHANGED, (long)1);

                            GDEBUG_STREAM("T1rho map, pmu time is " << t1pmap.headers_[offset].physiology_time_stamp[0]);
                        }
                    }
                }
            }

            Gadgetron::scal( (float)(scaling_factor_t1map.value()), t1map.data_);
            Gadgetron::scal( (float)(scaling_factor_t2map.value()), t2map.data_);
            if (this->has_t1p_mapping_) Gadgetron::scal((float)(scaling_factor_t1pmap.value()), t1pmap.data_);
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrParametricMappingGadget::fill_multi_map_header(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    // TODO: This should no longer be called, as it is replaced by perform_multi_mapping()
    int CmrParametricSashaT1T2MappingGadget::perform_mapping(IsmrmrdImageArray& data, IsmrmrdImageArray& map, IsmrmrdImageArray& para, IsmrmrdImageArray& map_sd, IsmrmrdImageArray& para_sd)
    {
        try
        {
            if (perform_timing.value()) { gt_timer_.start("CmrParametricSashaT1T2MappingGadget::perform_mapping"); }

            GDEBUG_CONDITION_STREAM(verbose.value(), "CmrParametricSashaT1T2MappingGadget::perform_mapping(...) starts ... ");

            size_t RO  = data.data_.get_size(0);
            size_t E1  = data.data_.get_size(1);
            size_t E2  = data.data_.get_size(2);
            size_t CHA = data.data_.get_size(3);
            size_t N   = data.data_.get_size(4);
            size_t S   = data.data_.get_size(5);
            size_t SLC = data.data_.get_size(6);

            size_t ro, e1, s, slc, p;

            GADGET_CHECK_RETURN(E2  == 1,                          GADGET_FAIL);
            GADGET_CHECK_RETURN(CHA == 1,                          GADGET_FAIL);
            GADGET_CHECK_RETURN(this->prep_times_ts_.size()  >= N, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->prep_times_t2p_.size() >= N, GADGET_FAIL);

            hoNDArray<float> mag;
            Gadgetron::abs(data.data_, mag);

            size_t curr_slc = data.headers_(0).slice;

            std::ostringstream str;
            str << "_slice_" << curr_slc;
            std::string slc_str = str.str();

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array(mag, debug_folder_full_path_ + "CmrParametricT1SRMapping_data_mag_" + slc_str);
            }

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<float> buf(N);
                memcpy(buf.begin(), &this->prep_times_ts_[0], sizeof(float)*N);
                gt_exporter_.export_array(buf, debug_folder_full_path_ + "CmrParametricT1SRMapping_prep_times_ts_" + slc_str);
            }

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<float> buf(N);
                memcpy(buf.begin(), &this->prep_times_t2p_[0], sizeof(float)*N);
                gt_exporter_.export_array(buf, debug_folder_full_path_ + "CmrParametricT1SRMapping_prep_times_t2p_" + slc_str);
            }

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<float> buf(this->t2p_rf_duration_.size());
                memcpy(&buf[0], &this->t2p_rf_duration_[0], buf.get_number_of_bytes());
                gt_exporter_.export_array(buf, debug_folder_full_path_ + "CmrParametricT1SRMapping_t2p_rf_duration_" + slc_str);
            }

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<float> buf(1);
                buf[0] = this->time_t2p_to_center_kspace_;
                gt_exporter_.export_array(buf, debug_folder_full_path_ + "CmrParametricT1SRMapping_time_t2p_to_center_kspace_" + slc_str);
            }

            bool need_sd_map = send_sd_map.value();

            Gadgetron::GadgetronTimer gt_timer(false);

            // -------------------------------------------------------------
            // set mapping parameters

            Gadgetron::CmrSashaT1T2Mapping<float> t1t2_sasha;

            t1t2_sasha.fill_holes_in_maps_ = perform_hole_filling.value();
            t1t2_sasha.max_size_of_holes_  = max_size_hole.value();
            t1t2_sasha.compute_SD_maps_    = need_sd_map;

            t1t2_sasha.ti_.resize(N*2+2, 0);
            memcpy(&(t1t2_sasha.ti_)[0],     &this->prep_times_ts_[0],          sizeof(float)*N);
            memcpy(&(t1t2_sasha.ti_)[N],     &this->prep_times_t2p_[0],         sizeof(float)*N);
            memcpy(&(t1t2_sasha.ti_)[N*2],   &this->t2p_rf_duration_,           sizeof(float)*1);
            memcpy(&(t1t2_sasha.ti_)[N*2+1], &this->time_t2p_to_center_kspace_, sizeof(float)*1);

            for (size_t n = 0; n < t1t2_sasha.ti_.size(); n++)
            {
                GDEBUG_STREAM("t1t2_sasha.ti[" << n << "] = " << t1t2_sasha.ti_[n]);
            }

            // set the anchor image TS
            size_t anchor_ind = this->anchor_image_index.value();
            if (anchor_ind < N)
            {
                t1t2_sasha.ti_[anchor_ind] = this->anchor_TS.value();
            }

            t1t2_sasha.data_.create(RO, E1, N, S, SLC, mag.begin());

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array(t1t2_sasha.data_, debug_folder_full_path_ + "CmrParametricT1SRMapping_t1t2_sasha_data_" + slc_str);
            }

            t1t2_sasha.max_iter_       = max_iter.value();
            t1t2_sasha.thres_fun_      = thres_func.value();
            t1t2_sasha.max_map_value_  = max_T1.value();

            t1t2_sasha.verbose_        = verbose.value();
            t1t2_sasha.debug_folder_   = debug_folder_full_path_;
            t1t2_sasha.perform_timing_ = perform_timing.value();

            // -------------------------------------------------------------
            // compute mask if needed
            if (mapping_with_masking.value())
            {
                t1t2_sasha.mask_for_mapping_.create(RO, E1, SLC);

                // get the image with longest TS time
                hoNDArray<float> mag_longest_TS;
                mag_longest_TS.create(RO, E1, SLC);

                for (slc = 0; slc < SLC; slc++)
                {
                    size_t ind = N - 1;
                    if (anchor_ind < N)
                    {
                        ind = anchor_ind;
                    }
                    else
                    {
                        float max_ts = this->prep_times_ts_[0];
                        for (size_t n = 1; n < this->prep_times_ts_.size(); n++)
                        {
                            if(this->prep_times_ts_[n]>max_ts)
                            {
                                max_ts = this->prep_times_ts_[n];
                                ind = n;
                            }
                        }
                    }

                    memcpy(&mag_longest_TS(0, 0, slc), &mag(0, 0, ind, 0, slc), sizeof(float)*RO*E1);
                }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(mag_longest_TS, debug_folder_full_path_ + "CmrParametricT1SRMapping_mag_longest_TS_" + slc_str);
                }

                double scale_factor = 1.0;
                if (data.meta_[0].length(GADGETRON_IMAGE_SCALE_RATIO) > 0)
                {
                    scale_factor = data.meta_[0].as_double(GADGETRON_IMAGE_SCALE_RATIO);
                }

                GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, find incoming image has scale factor of " << scale_factor);

                if (perform_timing.value()) { gt_timer.start("CmrParametricSashaT1T2MappingGadget::compute_mask_for_mapping"); }
                this->compute_mask_for_mapping(mag, t1t2_sasha.mask_for_mapping_, (float)scale_factor);
                if (perform_timing.value()) { gt_timer.stop(); }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(t1t2_sasha.mask_for_mapping_, debug_folder_full_path_ + "CmrParametricT1SRMapping_mask_for_mapping_" + slc_str);
                }
            }

            // -------------------------------------------------------------
            // perform mapping

            if (perform_timing.value()) { gt_timer.start("CmrParametricSashaT1T2MappingGadget, t1t2_sasha.perform_parametric_mapping"); }
            t1t2_sasha.perform_parametric_mapping();
            if (perform_timing.value()) { gt_timer.stop(); }

            size_t num_para = t1t2_sasha.get_num_of_paras();

            // -------------------------------------------------------------
            // get the results

            map.data_.create(RO, E1, E2, CHA, 1, S, SLC);
            Gadgetron::clear(map.data_);
            map.headers_.create(1, S, SLC);
            map.meta_.resize(S*SLC);

            para.data_.create(RO, E1, E2, CHA, num_para, S, SLC);
            Gadgetron::clear(para.data_);
            para.headers_.create(num_para, S, SLC);
            para.meta_.resize(num_para*S*SLC);

            if (need_sd_map)
            {
                map_sd.data_.create(RO, E1, E2, CHA, 1, S, SLC);
                Gadgetron::clear(map_sd.data_);
                map_sd.headers_.create(1, S, SLC);
                map_sd.meta_.resize(S*SLC);

                para_sd.data_.create(RO, E1, E2, CHA, num_para, S, SLC);
                Gadgetron::clear(para_sd.data_);
                para_sd.headers_.create(num_para, S, SLC);
                para_sd.meta_.resize(num_para*S*SLC);
            }

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (e1 = 0; e1 < E1; e1++)
                    {
                        for (ro = 0; ro < RO; ro++)
                        {
                            map.data_(ro, e1, 0, 0, 0, s, slc) = t1t2_sasha.map_(ro, e1, s, slc);

                            if (need_sd_map)
                            {
                                map_sd.data_(ro, e1, 0, 0, 0, s, slc) = t1t2_sasha.sd_map_(ro, e1, s, slc);
                            }

                            for (p = 0; p < num_para; p++)
                            {
                                para.data_(ro, e1, 0, 0, p, s, slc) = t1t2_sasha.para_(ro, e1, p, s, slc);

                                if (need_sd_map)
                                {
                                    para_sd.data_(ro, e1, 0, 0, p, s, slc) = t1t2_sasha.sd_para_(ro, e1, p, s, slc);
                                }
                            }
                        }
                    }

                    size_t slc_ind = data.headers_(0, s, slc).slice;

                    // Use the T1/T2 maps from para instead
                    map.headers_(0, s, slc) = data.headers_(0, s, slc);
                    map.headers_(0, s, slc).image_index = 1 + slc_ind;
                    map.headers_(0, s, slc).image_series_index = 11;
                    map.meta_[s+slc*S] = data.meta_[s + slc*S];
                    map.meta_[s + slc*S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1MAP);
                    map.meta_[s + slc*S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1MAP);
                    map.meta_[s + slc*S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1MAP);

                    if (need_sd_map)
                    {
                        map_sd.headers_(0, s, slc) = data.headers_(0, s, slc);
                        map_sd.headers_(0, s, slc).image_index = 1 + slc_ind;
                        map_sd.headers_(0, s, slc).image_series_index = 12;
                        map_sd.meta_[s + slc*S] = data.meta_[s + slc*S];
                        map_sd.meta_[s + slc*S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1SDMAP);
                        map_sd.meta_[s + slc*S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1SDMAP);
                        map_sd.meta_[s + slc*S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1SDMAP);
                    }

                    for (p = 0; p < num_para; p++)
                    {
                        para.headers_(p, s, slc) = data.headers_(0, s, slc);
                        para.headers_(p, s, slc).image_index = 1 + p + slc_ind*num_para;
                        para.meta_[p + s*num_para + slc*num_para*S] = data.meta_[s + slc*S];

                        if (need_sd_map)
                        {
                            para_sd.headers_(p, s, slc) = data.headers_(0, s, slc);
                            para_sd.headers_(p, s, slc).image_index = 1 + p + slc_ind*num_para;
                            para_sd.meta_[p + s*num_para + slc*num_para*S] = data.meta_[s + slc*S];
                        }
                    }
                }
            }

            // -------------------------------------------------------------

            if (perform_timing.value()) { gt_timer_.stop(); }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrParametricSashaT1T2MappingGadget::perform_mapping(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int CmrParametricSashaT1T2MappingGadget::perform_multi_mapping(IsmrmrdImageArray& data, IsmrmrdImageArray& t1map, IsmrmrdImageArray& t2map, IsmrmrdImageArray& para, IsmrmrdImageArray& map_sd, IsmrmrdImageArray& para_sd)
    {
        try
        {
            if (perform_timing.value()) { gt_timer_.start("CmrParametricSashaT1T2MappingGadget::perform_multi_mapping(t1, t2)"); }

            GDEBUG_CONDITION_STREAM(verbose.value(), "CmrParametricSashaT1T2MappingGadget::perform_multi_mapping(t1, t2) starts ... ");

            size_t RO  = data.data_.get_size(0);
            size_t E1  = data.data_.get_size(1);
            size_t E2  = data.data_.get_size(2);
            size_t CHA = data.data_.get_size(3);
            size_t N   = data.data_.get_size(4);
            size_t S   = data.data_.get_size(5);
            size_t SLC = data.data_.get_size(6);

            size_t ro, e1, s, slc, p;

            GADGET_CHECK_RETURN(E2  == 1,                          GADGET_FAIL);
            GADGET_CHECK_RETURN(CHA == 1,                          GADGET_FAIL);
            GADGET_CHECK_RETURN(this->prep_times_ts_.size()  >= N, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->prep_times_t2p_.size() >= N, GADGET_FAIL);

            hoNDArray<float> mag;
            Gadgetron::abs(data.data_, mag);

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array(mag, debug_folder_full_path_ + "CmrParametricT1SRMapping_data_mag");
            }

            bool need_sd_map = send_sd_map.value();

            Gadgetron::GadgetronTimer gt_timer(false);

            // ===========================================================
            // set mapping parameters

            Gadgetron::CmrSashaT1T2Mapping<float> t1t2_sasha;

            t1t2_sasha.fill_holes_in_maps_ = perform_hole_filling.value();
            t1t2_sasha.max_size_of_holes_  = max_size_hole.value();
            t1t2_sasha.compute_SD_maps_    = need_sd_map;

            // ti_ is a vector containing TS, TE, TSL, time-t2p-to-center, and T2p duration, of size 4N+1, where N is the size of y.
            // For the 'i'th measurement in y:
            //   x[i]     is the sat recovery time
            //   x[i+ N]  is the T2p time
            //   x[i+2N]  is the T1p time
            //   x[i+3N]  is the t2pRfDuration
            //   x[end]   is the timeT2pToCenter
            t1t2_sasha.ti_.resize(N*4 + 1, 0);
            memcpy(&(t1t2_sasha.ti_)[0],     &this->prep_times_ts_[  0],        sizeof(float) * N);
            memcpy(&(t1t2_sasha.ti_)[N],     &this->prep_times_t2p_[ 0],        sizeof(float) * N);
            memcpy(&(t1t2_sasha.ti_)[N*2],   &this->prep_times_t1p_[ 0],        sizeof(float) * N);
            memcpy(&(t1t2_sasha.ti_)[N*3],   &this->t2p_rf_duration_[0],        sizeof(float) * N);
            memcpy(&(t1t2_sasha.ti_)[N*4],   &this->time_t2p_to_center_kspace_, sizeof(float) * 1);

            GDEBUG_STREAM("======================================");
            for (size_t n = 0; n < this->prep_times_ts_.size(); n++)
            {
                GDEBUG_STREAM("this->prep_times_ts_[" << n << "] = " << this->prep_times_ts_[n]);
            }
            GDEBUG_STREAM("======================================");
            for (size_t n = 0; n < this->prep_times_t2p_.size(); n++)
            {
                GDEBUG_STREAM("this->prep_times_t2p_[" << n << "] = " << this->prep_times_t2p_[n]);
            }

            for (size_t n = 0; n < t1t2_sasha.ti_.size(); n++)
            {
                GDEBUG_STREAM("t1t2_sasha.ti[" << n << "] = " << t1t2_sasha.ti_[n]);
            }

            size_t curr_slc = data.headers_(0).slice;

            std::ostringstream str;
            str << "_slice_" << curr_slc;

            std::string slc_str = str.str();

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<float> buf(N);
                memcpy(buf.begin(), &this->prep_times_ts_[0], sizeof(float)*N);
                gt_exporter_.export_array(buf, debug_folder_full_path_ + "CmrParametricT1SRMapping_prep_times_ts_" + slc_str);
            }

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<float> buf(N);
                memcpy(buf.begin(), &this->prep_times_t2p_[0], sizeof(float)*N);
                gt_exporter_.export_array(buf, debug_folder_full_path_ + "CmrParametricT1SRMapping_prep_times_t2p_" + slc_str);
            }

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<float> buf(this->t2p_rf_duration_.size());
                memcpy(&buf[0], &this->t2p_rf_duration_[0], buf.get_number_of_bytes());
                gt_exporter_.export_array(buf, debug_folder_full_path_ + "CmrParametricT1SRMapping_t2p_rf_duration_" + slc_str);
            }

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<float> buf(1);
                buf[0] = this->time_t2p_to_center_kspace_;
                gt_exporter_.export_array(buf, debug_folder_full_path_ + "CmrParametricT1SRMapping_time_t2p_to_center_kspace_" + slc_str);
            }
          
            // set the anchor image TS
            size_t anchor_ind = this->anchor_image_index.value();
            if (anchor_ind < N)
            {
                t1t2_sasha.ti_[anchor_ind] = this->anchor_TS.value();
            }

            t1t2_sasha.data_.create(RO, E1, N, S, SLC, mag.begin());

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array(t1t2_sasha.data_, debug_folder_full_path_ + "CmrParametricT1SRMapping_t1t2_sasha_data_" + slc_str);
            }

            t1t2_sasha.max_iter_       = max_iter.value();
            t1t2_sasha.thres_fun_      = thres_func.value();
            t1t2_sasha.max_map_value_  = max_T1.value();

            t1t2_sasha.verbose_        = verbose.value();
            t1t2_sasha.debug_folder_   = debug_folder_full_path_;
            t1t2_sasha.perform_timing_ = perform_timing.value();

            // -------------------------------------------------------------
            // compute mask if needed
            if (mapping_with_masking.value())
            {
                t1t2_sasha.mask_for_mapping_.create(RO, E1, SLC);

                // get the image with longest TS time
                hoNDArray<float> mag_longest_TS;
                mag_longest_TS.create(RO, E1, SLC);

                for (slc = 0; slc < SLC; slc++)
                {
                    size_t ind = N - 1;
                    if (anchor_ind < N)
                    {
                        ind = anchor_ind;
                    }
                    else
                    {
                        float max_ts = this->prep_times_ts_[0];
                        for (size_t n = 1; n < this->prep_times_ts_.size(); n++)
                        {
                            if(this->prep_times_ts_[n]>max_ts)
                            {
                                max_ts = this->prep_times_ts_[n];
                                ind = n;
                            }
                        }
                    }

                    memcpy(&mag_longest_TS(0, 0, slc), &mag(0, 0, ind, 0, slc), sizeof(float)*RO*E1);
                }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(mag_longest_TS, debug_folder_full_path_ + "CmrParametricT1SRMapping_mag_longest_TS_" + slc_str);
                }

                double scale_factor = 1.0;
                if (data.meta_[0].length(GADGETRON_IMAGE_SCALE_RATIO) > 0)
                {
                    scale_factor = data.meta_[0].as_double(GADGETRON_IMAGE_SCALE_RATIO);
                }

                GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, find incoming image has scale factor of " << scale_factor);

                if (perform_timing.value()) { gt_timer.start("CmrParametricSashaT1T2MappingGadget::compute_mask_for_mapping"); }
                this->compute_mask_for_mapping(mag, t1t2_sasha.mask_for_mapping_, (float)scale_factor);
                if (perform_timing.value()) { gt_timer.stop(); }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(t1t2_sasha.mask_for_mapping_, debug_folder_full_path_ + "CmrParametricT1SRMapping_mask_for_mapping_" + slc_str);
                }
            }

            // -------------------------------------------------------------
            // perform mapping

            if (perform_timing.value()) { gt_timer.start("CmrParametricSashaT1T2MappingGadget, t1t2_sasha.perform_parametric_mapping"); }
            t1t2_sasha.perform_parametric_mapping();
            if (perform_timing.value()) { gt_timer.stop(); }

            size_t num_para = t1t2_sasha.get_num_of_paras();

            // -------------------------------------------------------------
            // get the results

            t1map.data_.create(RO, E1, E2, CHA, 1, S, SLC);
            Gadgetron::clear(t1map.data_);
            t1map.headers_.create(1, S, SLC);
            t1map.meta_.resize(S*SLC);

            t2map.data_.create(RO, E1, E2, CHA, 1, S, SLC);
            Gadgetron::clear(t2map.data_);
            t2map.headers_.create(1, S, SLC);
            t2map.meta_.resize(S*SLC);

            para.data_.create(RO, E1, E2, CHA, num_para, S, SLC);
            Gadgetron::clear(para.data_);
            para.headers_.create(num_para, S, SLC);
            para.meta_.resize(num_para*S*SLC);

            if (need_sd_map)
            {
                map_sd.data_.create(RO, E1, E2, CHA, 1, S, SLC);
                Gadgetron::clear(map_sd.data_);
                map_sd.headers_.create(1, S, SLC);
                map_sd.meta_.resize(S*SLC);

                para_sd.data_.create(RO, E1, E2, CHA, num_para, S, SLC);
                Gadgetron::clear(para_sd.data_);
                para_sd.headers_.create(num_para, S, SLC);
                para_sd.meta_.resize(num_para*S*SLC);
            }

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (e1 = 0; e1 < E1; e1++)
                    {
                        for (ro = 0; ro < RO; ro++)
                        {
                            // map.data_(ro, e1, 0, 0, 0, s, slc) = t1t2_sasha.map_(ro, e1, s, slc);

                            if (need_sd_map)
                            {
                                map_sd.data_(ro, e1, 0, 0, 0, s, slc) = t1t2_sasha.sd_map_(ro, e1, s, slc);
                            }

                            for (p = 0; p < num_para; p++)
                            {
                                para.data_(ro, e1, 0, 0, p, s, slc) = t1t2_sasha.para_(ro, e1, p, s, slc);

                                if (need_sd_map)
                                {
                                    para_sd.data_(ro, e1, 0, 0, p, s, slc) = t1t2_sasha.sd_para_(ro, e1, p, s, slc);
                                }
                            }
                            t1map.data_(ro, e1, 0, 0, 0, s, slc) = t1t2_sasha.para_(ro, e1, 1, s, slc);
                            t2map.data_(ro, e1, 0, 0, 0, s, slc) = t1t2_sasha.para_(ro, e1, 2, s, slc);
                        }
                    }

                    size_t slc_ind = data.headers_(0, s, slc).slice;

                    // Use the T1/T2 maps from para instead
                    t1map.headers_(0, s, slc) = data.headers_(0, s, slc);
                    t1map.headers_(0, s, slc).image_index = 1 + slc_ind;
                    t1map.headers_(0, s, slc).image_series_index = 11;
                    t1map.headers_(0, s, slc).repetition = 0;
                    t1map.headers_(0, s, slc).average = 0;
                    t1map.meta_[s+slc*S] = data.meta_[s + slc*S];
                    t1map.meta_[s + slc*S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1MAP);
                    t1map.meta_[s + slc*S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1MAP);
                    t1map.meta_[s + slc*S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1MAP);

                    // Use the T1/T2 maps from para instead
                    t2map.headers_(0, s, slc) = data.headers_(0, s, slc);
                    t2map.headers_(0, s, slc).image_index = 1 + slc_ind;
                    t2map.headers_(0, s, slc).image_series_index = 12;
                    t2map.headers_(0, s, slc).repetition = 0;
                    t2map.headers_(0, s, slc).average = 0;
                    t2map.meta_[s+slc*S] = data.meta_[s + slc*S];
                    t2map.meta_[s + slc*S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T2MAP);
                    t2map.meta_[s + slc*S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T2MAP);
                    t2map.meta_[s + slc*S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T2MAP);

                    if (need_sd_map)
                    {
                        map_sd.headers_(0, s, slc) = data.headers_(0, s, slc);
                        map_sd.headers_(0, s, slc).image_index = 1 + slc_ind;
                        map_sd.headers_(0, s, slc).image_series_index = 13;
                        map_sd.headers_(0, s, slc).repetition = 0;
                        map_sd.headers_(0, s, slc).average = 0;
                        map_sd.meta_[s + slc*S] = data.meta_[s + slc*S];
                        map_sd.meta_[s + slc*S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1SDMAP);
                        map_sd.meta_[s + slc*S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1SDMAP);
                        map_sd.meta_[s + slc*S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1SDMAP);
                    }

                    for (p = 0; p < num_para; p++)
                    {
                        para.headers_(p, s, slc) = data.headers_(0, s, slc);
                        para.headers_(p, s, slc).image_index = 1 + p + slc_ind*num_para;
                        para.headers_(p, s, slc).repetition = 0;
                        para.headers_(p, s, slc).average = 0;
                        para.meta_[p + s*num_para + slc*num_para*S] = data.meta_[s + slc*S];

                        if (need_sd_map)
                        {
                            para_sd.headers_(p, s, slc) = data.headers_(0, s, slc);
                            para_sd.headers_(p, s, slc).image_index = 1 + p + slc_ind*num_para;
                            para_sd.headers_(p, s, slc).repetition = 0;
                            para_sd.headers_(p, s, slc).average = 0;
                            para_sd.meta_[p + s*num_para + slc*num_para*S] = data.meta_[s + slc*S];
                        }
                    }
                }
            }

            // -------------------------------------------------------------

            if (perform_timing.value()) { gt_timer_.stop(); }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrParametricSashaT1T2MappingGadget::perform_multi_mapping(t1, t2) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int CmrParametricSashaT1T2MappingGadget::perform_multi_mapping(IsmrmrdImageArray& data, IsmrmrdImageArray& t1map, IsmrmrdImageArray& t2map, IsmrmrdImageArray& t1pmap, IsmrmrdImageArray& para, IsmrmrdImageArray& map_sd, IsmrmrdImageArray& para_sd)
    {
        try
        {
            if (perform_timing.value()) { gt_timer_.start("CmrParametricSashaT1T2MappingGadget::perform_multi_mapping(t1, t2, t1p)"); }

            GDEBUG_CONDITION_STREAM(verbose.value(), "CmrParametricSashaT1T2MappingGadget::perform_multi_mapping(t1, t2, t1p) starts ... ");

            size_t RO = data.data_.get_size(0);
            size_t E1 = data.data_.get_size(1);
            size_t E2 = data.data_.get_size(2);
            size_t CHA = data.data_.get_size(3);
            size_t N = data.data_.get_size(4);
            size_t S = data.data_.get_size(5);
            size_t SLC = data.data_.get_size(6);

            size_t ro, e1, s, slc, p;

            GADGET_CHECK_RETURN(E2 == 1, GADGET_FAIL);
            GADGET_CHECK_RETURN(CHA == 1, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->prep_times_ts_.size() >= N, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->prep_times_t2p_.size() >= N, GADGET_FAIL);
            GADGET_CHECK_RETURN(this->prep_times_t1p_.size() >= N, GADGET_FAIL);

            hoNDArray<float> mag;
            Gadgetron::abs(data.data_, mag);

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array(mag, debug_folder_full_path_ + "CmrParametricT1SRMapping_data_mag");
            }

            bool need_sd_map = send_sd_map.value();

            Gadgetron::GadgetronTimer gt_timer(false);

            // ===========================================================
            // set mapping parameters

            Gadgetron::CmrSashaT1T2T1pMapping<float> t1t2t1p_sasha;

            t1t2t1p_sasha.fill_holes_in_maps_ = perform_hole_filling.value();
            t1t2t1p_sasha.max_size_of_holes_ = max_size_hole.value();
            t1t2t1p_sasha.compute_SD_maps_ = need_sd_map;

            // ti_ is a vector containing TS, TE, TSL, time-t2p-to-center, and T2p duration, of size 4N+1, where N is the size of y.
            // For the 'i'th measurement in y:
            //   x[i]     is the sat recovery time
            //   x[i+ N]  is the T2p time
            //   x[i+2N]  is the T1p time
            //   x[i+3N]  is the t2pRfDuration
            //   x[end]   is the timeT2pToCenter
            t1t2t1p_sasha.ti_.resize(N*4 + 1, 0);
            memcpy(&(t1t2t1p_sasha.ti_)[0],     &this->prep_times_ts_[  0],        sizeof(float) * N);
            memcpy(&(t1t2t1p_sasha.ti_)[N],     &this->prep_times_t2p_[ 0],        sizeof(float) * N);
            memcpy(&(t1t2t1p_sasha.ti_)[N*2],   &this->prep_times_t1p_[ 0],        sizeof(float) * N);
            memcpy(&(t1t2t1p_sasha.ti_)[N*3],   &this->t2p_rf_duration_[0],        sizeof(float) * N);
            memcpy(&(t1t2t1p_sasha.ti_)[N*4],   &this->time_t2p_to_center_kspace_, sizeof(float) * 1);

            if (this->verbose.value())
            {
                for (size_t n = 0; n < t1t2t1p_sasha.ti_.size(); n++)
                {
                    GDEBUG_STREAM("t1t2t1p_sasha.ti[" << n << "] = " << t1t2t1p_sasha.ti_[n]);
                }
            }

            // set the anchor image TS
            size_t anchor_ind = this->anchor_image_index.value();
            if (anchor_ind < N)
            {
                t1t2t1p_sasha.ti_[anchor_ind] = this->anchor_TS.value();
            }

            t1t2t1p_sasha.data_.create(RO, E1, N, S, SLC, mag.begin());

            t1t2t1p_sasha.max_iter_ = max_iter.value();
            t1t2t1p_sasha.thres_fun_ = thres_func.value();
            t1t2t1p_sasha.max_map_value_ = max_T1.value();

            t1t2t1p_sasha.verbose_ = verbose.value();
            t1t2t1p_sasha.debug_folder_ = debug_folder_full_path_;
            t1t2t1p_sasha.perform_timing_ = perform_timing.value();

            // -------------------------------------------------------------
            // compute mask if needed
            if (mapping_with_masking.value())
            {
                t1t2t1p_sasha.mask_for_mapping_.create(RO, E1, SLC);

                // get the image with longest TS time
                hoNDArray<float> mag_longest_TS;
                mag_longest_TS.create(RO, E1, SLC);

                for (slc = 0; slc < SLC; slc++)
                {
                    size_t ind = N - 1;
                    if (anchor_ind < N)
                    {
                        ind = anchor_ind;
                    }
                    else
                    {
                        float max_ts = this->prep_times_ts_[0];
                        for (size_t n = 1; n < this->prep_times_ts_.size(); n++)
                        {
                            if (this->prep_times_ts_[n] > max_ts)
                            {
                                max_ts = this->prep_times_ts_[n];
                                ind = n;
                            }
                        }
                    }

                    memcpy(&mag_longest_TS(0, 0, slc), &mag(0, 0, ind, 0, slc), sizeof(float) * RO * E1);
                }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(mag_longest_TS, debug_folder_full_path_ + "CmrParametricT1SRMapping_mag_longest_TS");
                }

                double scale_factor = 1.0;
                if (data.meta_[0].length(GADGETRON_IMAGE_SCALE_RATIO) > 0)
                {
                    scale_factor = data.meta_[0].as_double(GADGETRON_IMAGE_SCALE_RATIO);
                }

                GDEBUG_STREAM("CmrParametricSashaT1T2MappingGadget, find incoming image has scale factor of " << scale_factor);

                if (perform_timing.value()) { gt_timer.start("CmrParametricSashaT1T2MappingGadget::compute_mask_for_mapping"); }
                this->compute_mask_for_mapping(mag, t1t2t1p_sasha.mask_for_mapping_, (float)scale_factor);
                if (perform_timing.value()) { gt_timer.stop(); }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(t1t2t1p_sasha.mask_for_mapping_, debug_folder_full_path_ + "CmrParametricT1SRMapping_mask_for_mapping");
                }
            }

            // -------------------------------------------------------------
            // perform mapping

            if (perform_timing.value()) { gt_timer.start("CmrParametricSashaT1T2MappingGadget, t1t2t1p_sasha.perform_parametric_mapping"); }
            t1t2t1p_sasha.perform_parametric_mapping();
            if (perform_timing.value()) { gt_timer.stop(); }

            size_t num_para = t1t2t1p_sasha.get_num_of_paras();

            // -------------------------------------------------------------
            // get the results

            t1map.data_.create(RO, E1, E2, CHA, 1, S, SLC);
            Gadgetron::clear(t1map.data_);
            t1map.headers_.create(1, S, SLC);
            t1map.meta_.resize(S * SLC);

            t2map.data_.create(RO, E1, E2, CHA, 1, S, SLC);
            Gadgetron::clear(t2map.data_);
            t2map.headers_.create(1, S, SLC);
            t2map.meta_.resize(S * SLC);

            t1pmap.data_.create(RO, E1, E2, CHA, 1, S, SLC);
            Gadgetron::clear(t1pmap.data_);
            t1pmap.headers_.create(1, S, SLC);
            t1pmap.meta_.resize(S * SLC);

            para.data_.create(RO, E1, E2, CHA, num_para, S, SLC);
            Gadgetron::clear(para.data_);
            para.headers_.create(num_para, S, SLC);
            para.meta_.resize(num_para * S * SLC);

            if (need_sd_map)
            {
                map_sd.data_.create(RO, E1, E2, CHA, 1, S, SLC);
                Gadgetron::clear(map_sd.data_);
                map_sd.headers_.create(1, S, SLC);
                map_sd.meta_.resize(S * SLC);

                para_sd.data_.create(RO, E1, E2, CHA, num_para, S, SLC);
                Gadgetron::clear(para_sd.data_);
                para_sd.headers_.create(num_para, S, SLC);
                para_sd.meta_.resize(num_para * S * SLC);
            }

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (e1 = 0; e1 < E1; e1++)
                    {
                        for (ro = 0; ro < RO; ro++)
                        {
                            // map.data_(ro, e1, 0, 0, 0, s, slc) = t1t2t1p_sasha.map_(ro, e1, s, slc);

                            if (need_sd_map)
                            {
                                map_sd.data_(ro, e1, 0, 0, 0, s, slc) = t1t2t1p_sasha.sd_map_(ro, e1, s, slc);
                            }

                            for (p = 0; p < num_para; p++)
                            {
                                para.data_(ro, e1, 0, 0, p, s, slc) = t1t2t1p_sasha.para_(ro, e1, p, s, slc);

                                if (need_sd_map)
                                {
                                    para_sd.data_(ro, e1, 0, 0, p, s, slc) = t1t2t1p_sasha.sd_para_(ro, e1, p, s, slc);
                                }
                            }
                            t1map.data_(ro, e1, 0, 0, 0, s, slc) = t1t2t1p_sasha.para_(ro, e1, 1, s, slc);
                            t2map.data_(ro, e1, 0, 0, 0, s, slc) = t1t2t1p_sasha.para_(ro, e1, 2, s, slc);
                            t1pmap.data_(ro, e1, 0, 0, 0, s, slc) = t1t2t1p_sasha.para_(ro, e1, 3, s, slc);
                        }
                    }

                    size_t slc_ind = data.headers_(0, s, slc).slice;

                    // Use the T1/T2 maps from para instead
                    t1map.headers_(0, s, slc) = data.headers_(0, s, slc);
                    t1map.headers_(0, s, slc).image_index = 1 + slc_ind;
                    t1map.headers_(0, s, slc).image_series_index = 11;
                    t1map.headers_(0, s, slc).repetition = 0;
                    t1map.headers_(0, s, slc).average = 0;
                    t1map.meta_[s + slc * S] = data.meta_[s + slc * S];
                    t1map.meta_[s + slc * S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1MAP);
                    t1map.meta_[s + slc * S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1MAP);
                    t1map.meta_[s + slc * S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1MAP);

                    // Use the T1/T2 maps from para instead
                    t2map.headers_(0, s, slc) = data.headers_(0, s, slc);
                    t2map.headers_(0, s, slc).image_index = 1 + slc_ind;
                    t2map.headers_(0, s, slc).image_series_index = 12;
                    t2map.headers_(0, s, slc).repetition = 0;
                    t2map.headers_(0, s, slc).average = 0;
                    t2map.meta_[s + slc * S] = data.meta_[s + slc * S];
                    t2map.meta_[s + slc * S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T2MAP);
                    t2map.meta_[s + slc * S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T2MAP);
                    t2map.meta_[s + slc * S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T2MAP);

                    t1pmap.headers_(0, s, slc) = data.headers_(0, s, slc);
                    t1pmap.headers_(0, s, slc).image_index = 1 + slc_ind;
                    t1pmap.headers_(0, s, slc).image_series_index = 14;
                    t1pmap.headers_(0, s, slc).repetition = 0;
                    t1pmap.headers_(0, s, slc).average = 0;
                    t1pmap.meta_[s + slc * S] = data.meta_[s + slc * S];
                    t1pmap.meta_[s + slc * S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1RHOMAP);
                    t1pmap.meta_[s + slc * S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1RHOMAP);
                    t1pmap.meta_[s + slc * S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1RHOMAP);

                    if (need_sd_map)
                    {
                        map_sd.headers_(0, s, slc) = data.headers_(0, s, slc);
                        map_sd.headers_(0, s, slc).image_index = 1 + slc_ind;
                        map_sd.headers_(0, s, slc).image_series_index = 13;
                        map_sd.headers_(0, s, slc).repetition = 0;
                        map_sd.headers_(0, s, slc).average = 0;
                        map_sd.meta_[s + slc * S] = data.meta_[s + slc * S];
                        map_sd.meta_[s + slc * S].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_T1SDMAP);
                        map_sd.meta_[s + slc * S].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_T1SDMAP);
                        map_sd.meta_[s + slc * S].append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_IMAGE_T1SDMAP);
                    }

                    for (p = 0; p < num_para; p++)
                    {
                        para.headers_(p, s, slc) = data.headers_(0, s, slc);
                        para.headers_(p, s, slc).image_index = 1 + p + slc_ind * num_para;
                        para.headers_(p, s, slc).repetition = 0;
                        para.headers_(p, s, slc).average = 0;
                        para.meta_[p + s * num_para + slc * num_para * S] = data.meta_[s + slc * S];

                        if (need_sd_map)
                        {
                            para_sd.headers_(p, s, slc) = data.headers_(0, s, slc);
                            para_sd.headers_(p, s, slc).image_index = 1 + p + slc_ind * num_para;
                            para_sd.headers_(p, s, slc).repetition = 0;
                            para_sd.headers_(p, s, slc).average = 0;
                            para_sd.meta_[p + s * num_para + slc * num_para * S] = data.meta_[s + slc * S];
                        }
                    }
                }
            }

            // -------------------------------------------------------------

            if (perform_timing.value()) { gt_timer_.stop(); }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrParametricSashaT1T2MappingGadget::perform_multi_mapping(t1, t2, t1p) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int CmrParametricSashaT1T2MappingGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "CmrParametricSashaT1T2MappingGadget - close(flags) : " << flags);

        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (flags != 0)
        {
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(CmrParametricSashaT1T2MappingGadget)

}
