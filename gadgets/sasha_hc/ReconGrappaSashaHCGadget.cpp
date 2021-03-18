#include "ReconGrappaSashaHCGadget.h"
#include <mri_core_grappa.h>

namespace Gadgetron {

    ReconGrappaSashaHCGadget::ReconGrappaSashaHCGadget() : GenericReconCartesianGrappaGadget()
    {
    }

    ReconGrappaSashaHCGadget::~ReconGrappaSashaHCGadget()
    {
    }

    int ReconGrappaSashaHCGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        // Initialize storage of GRAPPA-ed k-space data from previous encoding space
        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        std::vector<size_t> dim;
        recon_bit_->rbit_[0].data_.data_.get_dimensions(dim);
        unaliased_uncombined_enc0_.create(dim);

        // Call base class
        return ( GenericReconCartesianGrappaGadget::process(m1) );
    }

    void ReconGrappaSashaHCGadget::perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        // Kelvin: This is largely copied from GenericReconCartesianGrappaGadget, with additions marked with "// Kelvin"
        try
        {
            typedef std::complex<float> T;

            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);
            size_t dstCHA = recon_bit.data_.data_.get_size(3);
            size_t N = recon_bit.data_.data_.get_size(4);
            size_t S = recon_bit.data_.data_.get_size(5);
            size_t SLC = recon_bit.data_.data_.get_size(6);

            hoNDArray< std::complex<float> >& src = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            size_t unmixingCoeff_CHA = recon_obj.unmixing_coeff_.get_size(3);

            size_t convkRO = recon_obj.kernel_.get_size(0);
            size_t convkE1 = recon_obj.kernel_.get_size(1);
            size_t convkE2 = recon_obj.kernel_.get_size(2);

            recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << e;
                std::string suffix = os.str();
                gt_exporter_.export_array_complex(recon_bit.data_.data_, debug_folder_full_path_ + "data_src_" + suffix);
            }

            // compute aliased images
            data_recon_buf_.create(RO, E1, E2, dstCHA, N, S, SLC);

            if (E2>1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }

            // SNR unit scaling
            float effective_acce_factor(1), snr_scaling_ratio(1);
            this->compute_snr_scaling_factor(recon_bit, effective_acce_factor, snr_scaling_ratio);
            if (effective_acce_factor > 1)
            {
                // since the grappa in gadgetron is doing signal preserving scaling, to perserve noise level, we need this compensation factor
                double grappaKernelCompensationFactor = 1.0 / (acceFactorE1_[e] * acceFactorE2_[e]);
                Gadgetron::scal((float)(grappaKernelCompensationFactor*snr_scaling_ratio), complex_im_recon_buf_);

                if (this->verbose.value()) GDEBUG_STREAM("GenericReconCartesianGrappaGadget, grappaKernelCompensationFactor*snr_scaling_ratio : " << grappaKernelCompensationFactor*snr_scaling_ratio);
            }

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << e;
                std::string suffix = os.str();
                gt_exporter_.export_array_complex(complex_im_recon_buf_, debug_folder_full_path_ + "aliasedIm_" + suffix);
            }

            // unwrapping

            long long num = N*S*SLC;

            long long ii;

#pragma omp parallel default(none) private(ii) shared(num, N, S, RO, E1, E2, srcCHA, convkRO, convkE1, convkE2, ref_N, ref_S, recon_obj, recon_bit, dstCHA, unmixingCoeff_CHA, e) if(num>1)
            {
#pragma omp for 
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (N*S);
                    size_t s = (ii - slc*N*S) / N;
                    size_t n = ii - slc*N*S - s*N;

                    // combined channels
                    T* pIm = &(complex_im_recon_buf_(0, 0, 0, 0, n, s, slc));

                    size_t usedN = n;
                    if (n >= ref_N) usedN = ref_N - 1;

                    size_t usedS = s;
                    if (s >= ref_S) usedS = ref_S - 1;

                    T* pUnmix = &(recon_obj.unmixing_coeff_(0, 0, 0, 0, usedN, usedS, slc));

                    T* pRes = &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc));
                    hoNDArray< std::complex<float> > res(RO, E1, E2, 1, pRes);

                    hoNDArray< std::complex<float> > unmixing(RO, E1, E2, unmixingCoeff_CHA, pUnmix);
                    hoNDArray< std::complex<float> > aliasedIm(RO, E1, E2, ((unmixingCoeff_CHA<=srcCHA) ? unmixingCoeff_CHA : srcCHA), 1, pIm);

                    // -------------------------------------
                    // Kelvin: START of SASHA-HC additions
                    // -------------------------------------
                    // Standard logic.  Don't do this yet, as we're doing different logic for first/subsequent encoding spaces
                    // Gadgetron::apply_unmix_coeff_aliased_image_3D(aliasedIm, unmixing, res);

                    // Determine the min/max acquired phase encode line
                    size_t minE1 = 65536;
                    size_t maxE1 = 0;

                    for (size_t e1 = 0; e1 < E1; e1++)
                    {
                        if (std::abs(recon_bit.data_.data_(RO / 2, e1, 0, 0, n, s)) > 0)
                        {
                            if (e1 < minE1)
                            {
                                minE1 = e1;
                            }

                            if (e1 > maxE1)
                            {
                                maxE1 = e1;
                            }
                        }
                    }
                    GDEBUG_STREAM("e: " << e << " minE1: " << minE1 << " maxE1: " << maxE1);

                    if (e == 0)
                    {
                        // ----------------------------------------------------------------------------
                        // For the first encoding space, just store uncombined (unaliased) k-space data
                        // ----------------------------------------------------------------------------
                        // Do normal unaliasing first
                        Gadgetron::apply_unmix_coeff_aliased_image_3D(aliasedIm, unmixing, res);

                        hoNDArray< std::complex<float> > unaliased_enc0(RO, E1, E2, unmixingCoeff_CHA, &(unaliased_uncombined_enc0_(0, 0, 0, 0, n, s, slc)) );

                        // Image-space GRAPPA and iFFT directly into unaliased matrix
                        Gadgetron::multiply(aliasedIm, unmixing, unaliased_enc0);
                        Gadgetron::hoNDFFT<float>::instance()->fft2c(unaliased_enc0);
                    } else {
                        // ----------------------------------------------------------------------------
                        // For subsequent encoding spaces, copy specific k-space data from the first
                        // ----------------------------------------------------------------------------
                        // Reference to data from first encoding space
                        hoNDArray< std::complex<float> > unaliased_enc0(RO, E1, E2, unmixingCoeff_CHA, &(unaliased_uncombined_enc0_(0, 0, 0, 0, n, s, slc)) );
                        hoNDArray< std::complex<float> > combined(      RO, E1, E2, 1,                 &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc)) );

                        // Debug export
                        if (!debug_folder_full_path_.empty())
                        {
                            std::stringstream os;
                            os << "enc_" << e;
                            std::string suffix = os.str();
                            gt_exporter_.export_array_complex(unaliased_enc0, debug_folder_full_path_ + "buffer0_" + suffix);
                        }

                        // Create a buffer object of the same size
                        //std::vector<size_t> dim;
                        //aliasedIm.get_dimensions(dim);
                        hoNDArray<T> buffer;
                        buffer.create(RO, E1, E2, unmixingCoeff_CHA);

                        // Image-space GRAPPA and iFFT into a buffer
                        Gadgetron::multiply(aliasedIm, unmixing, buffer);
                        Gadgetron::hoNDFFT<float>::instance()->fft2c(buffer);

                        // Debug export
                        if (!debug_folder_full_path_.empty())
                        {
                            std::stringstream os;
                            os << "enc_" << e;
                            std::string suffix = os.str();
                            gt_exporter_.export_array_complex(buffer, debug_folder_full_path_ + "bufferBefore_" + suffix);
                        }

                        // Copy relevant lines to from first space
                        size_t ro, e1, e2, cha;
                        for (cha = 0; cha < unmixingCoeff_CHA; cha++)
                        {
                            for (e2 = 0; e2 < E2; e2++)
                            {
                                for (e1 = maxE1; e1 < E1; e1++)
                                {
                                    for (ro = 0; ro < RO; ro++)
                                    {
                                        buffer(ro, e1, e2, cha) = unaliased_enc0(ro, e1, e2, cha);
                                    }
                                }
                            }
                        }

                        // Debug export
                        if (!debug_folder_full_path_.empty())
                        {
                            std::stringstream os;
                            os << "enc_" << e;
                            std::string suffix = os.str();
                            gt_exporter_.export_array_complex(buffer, debug_folder_full_path_ + "bufferAfter_" + suffix);
                        }

                        // FFT back and store in the actual result
                        Gadgetron::hoNDFFT<float>::instance()->ifft2c(buffer);

                        // Debug export
                        if (!debug_folder_full_path_.empty())
                        {
                            std::stringstream os;
                            os << "enc_" << e;
                            std::string suffix = os.str();
                            gt_exporter_.export_array_complex(buffer, debug_folder_full_path_ + "bufferAfterFft_" + suffix);
                        }

                        Gadgetron::sum_over_dimension(buffer, res, 3);
                    } // else if (e == 0)

                    // Debug export
                    if (!debug_folder_full_path_.empty())
                    {
                        std::stringstream os;
                        os << "enc_" << e << "_n_" << n << "_s_" << s;
                        std::string suffix = os.str();
                        gt_exporter_.export_array_complex(res, debug_folder_full_path_ + "resAfter_" + suffix);
                    }

                    // -------------------------------------
                    // Kelvin: END of SASHA-HC additions
                    // -------------------------------------
                } // for (ii = 0; ii < num; ii++)
            }

            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << e;
                std::string suffix = os.str();
                gt_exporter_.export_array_complex(recon_obj.recon_res_.data_, debug_folder_full_path_ + "unwrappedIm_" + suffix);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconCartesianGrappaGadget::perform_unwrapping(...) ... ");
        }
    }

    GADGET_FACTORY_DECLARE( ReconGrappaSashaHCGadget )
}
