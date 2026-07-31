#include "DicomOrientationGadget.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "mri_core_utility.h"

namespace Gadgetron
{
    void DicomOrientationGadget::find_max_in_vec(float vector[3], float& elem, int& orientation)
    {
        if (std::abs(vector[2]) > std::abs(vector[1]))
        {
            elem = vector[2];
            orientation = 2;
            if (std::abs(vector[0]) > std::abs(vector[2]))
            {
                elem = vector[0];
                orientation = 0;
            }
        }
        else
        {
            elem = vector[1];
            orientation = 1;
            if (std::abs(vector[0]) > std::abs(vector[1]))
            {
                elem = vector[0];
                orientation = 0;
            }
        }
    }

    void DicomOrientationGadget::find_min_in_vec(float vector[3], float& elem, int& orientation)
    {
        if (std::abs(vector[2]) < std::abs(vector[1]))
        {
            elem = vector[2];
            orientation = 2;
            if (std::abs(vector[0]) < std::abs(vector[2]))
            {
                elem = vector[0];
                orientation = 0;
            }
        }
        else
        {
            elem = vector[1];
            orientation = 1;
            if (std::abs(vector[0]) < std::abs(vector[1]))
            {
                elem = vector[0];
                orientation = 0;
            }
        }
    }

    void DicomOrientationGadget::compute_norm_form_row_col(float rowvec[3], float colvec[3], float normvec[3])
    {
        normvec[0] = rowvec[1] * colvec[2] - rowvec[2] * colvec[1];
        normvec[1] = rowvec[2] * colvec[0] - rowvec[0] * colvec[2];
        normvec[2] = rowvec[0] * colvec[1] - rowvec[1] * colvec[0];
    }

    void DicomOrientationGadget::extract_main_orientation(float rowvec[3], float colvec[3], float normvec[3], float& value_norm, float& value_row, float& value_column, int& orientation_main, int& orientation_row)
    {
        float maxNorm(0); int maxOriNorm(0);
        find_max_in_vec(normvec, maxNorm, maxOriNorm);

        orientation_main = maxOriNorm;
        value_norm = maxNorm;

        float maxRow(0); int maxOriRow(0);
        float minRow(0); int minOriRow(0);
        find_max_in_vec(rowvec, maxRow, maxOriRow);
        find_min_in_vec(rowvec, minRow, minOriRow);

        float maxCol(0); int maxOriCol(0);
        float minCol(0); int minOriCol(0);
        find_max_in_vec(colvec, maxCol, maxOriCol);
        find_min_in_vec(colvec, minCol, minOriCol);

        float d_row(0), d_col(0);

        if (minOriRow == 2)
            d_row = std::abs(std::abs(rowvec[0]) - std::abs(rowvec[1]));
        else if (minOriRow == 1)
            d_row = std::abs(std::abs(rowvec[0]) - std::abs(rowvec[2]));
        else
            d_row = std::abs(std::abs(rowvec[2]) - std::abs(rowvec[1]));

        if (minOriCol == 2)
            d_col = std::abs(std::abs(colvec[0]) - std::abs(colvec[1]));
        else if (minOriCol == 1)
            d_col = std::abs(std::abs(colvec[0]) - std::abs(colvec[2]));
        else
            d_col = std::abs(std::abs(colvec[2]) - std::abs(colvec[1]));

        bool row_t = true, row_s = true, row_c = true;
        bool col_t = true, col_s = true, col_c = true;
        if (maxOriNorm == 2)
        {
            row_t = false;
            col_t = false;
        }
        else if (maxOriNorm == 1)
        {
            row_c = false;
            col_c = false;
        }
        else
        {
            row_s = false;
            col_s = false;
        }

        float oriRow(0), oriCol(0);
        int mainOriRow(0);

        if (d_row > d_col)
        {
            if (maxOriRow == 1)
            {
                col_c = false;
                oriRow = rowvec[1];
                mainOriRow = 1;
            }
            else if (maxOriRow == 0)
            {
                col_s = false;
                oriRow = rowvec[0];
                mainOriRow = 0;
            }
            else
            {
                col_t = false;
                oriRow = rowvec[2];
                mainOriRow = 2;
            }
        }
        else
        {
            if (maxOriCol == 1)
            {
                row_c = false;
                oriCol = colvec[1];
            }
            else if (maxOriCol == 0)
            {
                row_s = false;
                oriCol = colvec[0];
            }
            else
            {
                row_t = false;
                oriCol = colvec[2];
            }
        }

        if (d_row > d_col)
        {
            if (col_s) oriCol = colvec[0];
            if (col_c) oriCol = colvec[1];
            if (col_t) oriCol = colvec[2];
        }
        else
        {
            if (row_s) { oriRow = rowvec[0]; mainOriRow = 0; }
            if (row_c) { oriRow = rowvec[1]; mainOriRow = 1; }
            if (row_t) { oriRow = rowvec[2]; mainOriRow = 2; }
        }

        orientation_row = mainOriRow;
        value_row = oriRow;
        value_column = oriCol;
    }

    void DicomOrientationGadget::ismrmd_image_2d_to_dicom_compute_mirror_rotate(float* read_dir, float* phase_dir, float* slice_dir, float& rotate_angle, int& mirror)
    {
        float* colvec = read_dir;
        float* rowvec = phase_dir;
        float* normvec = slice_dir;

        float value_norm(0);
        float value_row(0);
        float value_column(0);
        int orientation_main(0);
        int orientation_row(0);

        extract_main_orientation(rowvec, colvec, normvec, value_norm, value_row, value_column, orientation_main, orientation_row);

        float mainElementRow = value_row;
        float mainElementColumn = value_column;
        int mainOrientation = orientation_main;
        int rowOrientation = orientation_row;

        bool isOrientedRow = false;
        bool isOrientedColumn = false;

        int desiredOrientationRow(0);

        switch (mainOrientation)
        {
        case 0:
            desiredOrientationRow = 1;
            break;
        case 1:
            desiredOrientationRow = 0;
            break;
        case 2:
            desiredOrientationRow = 0;
            break;
        }

        switch (mainOrientation)
        {
        case 0:
        case 1:
            if (rowOrientation == desiredOrientationRow)
            {
                isOrientedRow = mainElementRow > 0;
                isOrientedColumn = mainElementColumn < 0;
            }
            else
            {
                isOrientedRow = mainElementRow < 0;
                isOrientedColumn = mainElementColumn > 0;
            }
            break;
        case 2:
            isOrientedRow = mainElementRow > 0;
            isOrientedColumn = mainElementColumn > 0;
            break;
        }

        rotate_angle = 0;
        mirror = -1;
        if (rowOrientation == desiredOrientationRow)
        {
            if (isOrientedRow && !isOrientedColumn)
            {
                mirror = 1;
            }

            if (!isOrientedRow && isOrientedColumn)
            {
                mirror = 0;
            }

            if (!isOrientedRow && !isOrientedColumn)
            {
                rotate_angle = 180;
            }
        }
        else
        {
            if (mainOrientation == 0)
            {
                isOrientedRow = !isOrientedRow;
                isOrientedColumn = !isOrientedColumn;
            }

            compute_norm_form_row_col(rowvec, colvec, normvec);

            float maxNorm;
            int maxOriNorm;
            find_max_in_vec(normvec, maxNorm, maxOriNorm);

            if (maxNorm > 0)
            {
                isOrientedRow = !isOrientedRow;
                isOrientedColumn = !isOrientedColumn;
            }

            if (isOrientedRow && isOrientedColumn)
            {
                rotate_angle = 90;
                mirror = 1;
            }

            if (isOrientedRow && !isOrientedColumn)
            {
                rotate_angle = 90;
            }

            if (!isOrientedRow && isOrientedColumn)
            {
                rotate_angle = 270;
            }

            if (!isOrientedRow && !isOrientedColumn)
            {
                rotate_angle = 270;
                mirror = 1;
            }
        }
    }

    DicomOrientationGadget::DicomOrientationGadget(): BaseClass()
    {
    }

    DicomOrientationGadget::~DicomOrientationGadget()
    {
    }

    int DicomOrientationGadget::process_config(ACE_Message_Block* mb)
    {
        if (!debug_folder.value().empty())
        {
            Gadgetron::get_debug_folder_path(debug_folder.value(), debug_folder_full_path_);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is " << debug_folder_full_path_);

            boost::filesystem::path boost_folder_path(debug_folder_full_path_);
            try
            {
                boost::filesystem::create_directories(boost_folder_path);
            }
            catch (...)
            {
                GERROR("Error creating the debug folder.\n");
                return GADGET_FAIL;
            }
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is not set ... ");
        }

        ISMRMRD::IsmrmrdHeader h;
        deserialize(mb->rd_ptr(), h);
        xml = h;

        // get the vendor
        if (h.acquisitionSystemInformation)
        {
            if (h.acquisitionSystemInformation->systemVendor)
            {
                vendor_ = *h.acquisitionSystemInformation->systemVendor;
            }
        }

        if (h.measurementInformation)
        {
            if (h.measurementInformation->protocolName)
            {
                protocol_ = *h.measurementInformation->protocolName;
                boost::algorithm::to_lower(protocol_);
            }
        }

        return GADGET_OK;
    }

    int DicomOrientationGadget::process(GadgetContainerMessage<IsmrmrdImageArray>* m1)
    {
        auto& imagearr = *m1->getObjectPtr();

        const size_t X = imagearr.data_.get_size(0);
        const size_t Y = imagearr.data_.get_size(1);
        const size_t Z = imagearr.data_.get_size(2);
        const size_t CHA = imagearr.data_.get_size(3);
        const size_t N = imagearr.data_.get_size(4);
        const size_t S = imagearr.data_.get_size(5);
        const size_t LOC = imagearr.data_.get_size(6);
        const size_t expected_headers = N * S * LOC;

        if (verbose.value())
        {
            GDEBUG_STREAM("DicomOrientationGadget::process(image array) data dims [X Y Z CHA N S LOC] = ["
                          << X << " " << Y << " " << Z << " " << CHA << " " << N << " " << S << " " << LOC << "]");
            GDEBUG_STREAM("DicomOrientationGadget::process(image array) headers dims [N S LOC] = ["
                          << imagearr.headers_.get_size(0) << " " << imagearr.headers_.get_size(1) << " "
                          << imagearr.headers_.get_size(2) << "], meta size = " << imagearr.meta_.size()
                          << ", expected headers = " << expected_headers);
        }

        if (imagearr.headers_.get_number_of_elements() == 0)
        {
            GERROR_STREAM("DicomOrientationGadget::process, image array has no headers");
            return GADGET_FAIL;
        }

        if (imagearr.headers_.get_size(0) != N || imagearr.headers_.get_size(1) != S || imagearr.headers_.get_size(2) != LOC)
        {
            GERROR_STREAM("DicomOrientationGadget::process, image array header dimensions ["
                          << imagearr.headers_.get_size(0) << " " << imagearr.headers_.get_size(1) << " "
                          << imagearr.headers_.get_size(2) << "] do not match data dimensions ["
                          << N << " " << S << " " << LOC << "]");
            return GADGET_FAIL;
        }

        hoNDArray<std::complex<float>> output_data;
        bool output_initialized = false;

        for (size_t loc = 0; loc < LOC; loc++)
        {
            for (size_t s = 0; s < S; s++)
            {
                for (size_t n = 0; n < N; n++)
                {
                    const size_t meta_index = loc * N * S + s * N + n;
                    const bool has_meta = meta_index < imagearr.meta_.size();

                    if (verbose.value())
                    {
                        GDEBUG_STREAM("DicomOrientationGadget::process(image array) element [n s loc] = ["
                                      << n << " " << s << " " << loc << "], meta_index = " << meta_index
                                      << ", has_meta = " << has_meta);
                    }

                    auto image_header = imagearr.headers_(n, s, loc);
                    std::vector<size_t> single_image_dims{ X, Y, Z, CHA };
                    auto image_data = hoNDArray<std::complex<float>>(single_image_dims);
                    std::memcpy(image_data.get_data_ptr(), &imagearr.data_(0, 0, 0, 0, n, s, loc), X * Y * Z * CHA * sizeof(std::complex<float>));

                    if (verbose.value())
                    {
                        GDEBUG_STREAM("DicomOrientationGadget::process(image array) element dims before adjust [X Y Z CHA] = ["
                                      << image_data.get_size(0) << " " << image_data.get_size(1) << " "
                                      << image_data.get_size(2) << " " << image_data.get_size(3) << "]");
                    }

                    ISMRMRD::MetaContainer meta_container;
                    ISMRMRD::MetaContainer* meta_ptr = NULL;
                    if (has_meta)
                    {
                        meta_container = imagearr.meta_[meta_index];
                        meta_ptr = &meta_container;
                    }

                    bool should_adjust = !check_correct_image_orientation_flag.value();
                    if (check_correct_image_orientation_flag.value())
                    {
                        should_adjust = has_meta && meta_container.length(GADGETRON_CORRECT_IMAGE_ORIENTATION) > 0;
                    }

                    if (should_adjust)
                    {
                        if (verbose.value())
                        {
                            GDEBUG_STREAM("DicomOrientationGadget::process(image array) adjusting element [n s loc] = ["
                                          << n << " " << s << " " << loc << "]");
                        }

                        if (this->adjust_image_core(image_header, image_data, meta_ptr) != GADGET_OK)
                        {
                            GERROR_STREAM("DicomOrientationGadget::adjust_image_core failed for image array element [n s loc] = [" << n << " " << s << " " << loc << "]");
                            return GADGET_FAIL;
                        }

                        if (verbose.value())
                        {
                            GDEBUG_STREAM("DicomOrientationGadget::process(image array) element dims after adjust [X Y Z CHA] = ["
                                          << image_data.get_size(0) << " " << image_data.get_size(1) << " "
                                          << image_data.get_size(2) << " " << image_data.get_size(3) << "]");
                        }
                    }
                    else if (verbose.value())
                    {
                        GDEBUG_STREAM("DicomOrientationGadget::process(image array) skipping adjustment for element [n s loc] = ["
                                      << n << " " << s << " " << loc << "]");
                    }

                    if (!output_initialized)
                    {
                        output_data.create(image_data.get_size(0), image_data.get_size(1), image_data.get_size(2), image_data.get_size(3), N, S, LOC);
                        output_initialized = true;

                        if (verbose.value())
                        {
                            GDEBUG_STREAM("DicomOrientationGadget::process(image array) output buffer created with dims [X Y Z CHA N S LOC] = ["
                                          << output_data.get_size(0) << " " << output_data.get_size(1) << " "
                                          << output_data.get_size(2) << " " << output_data.get_size(3) << " "
                                          << output_data.get_size(4) << " " << output_data.get_size(5) << " "
                                          << output_data.get_size(6) << "]");
                        }
                    }
                    else if (output_data.get_size(0) != image_data.get_size(0) || output_data.get_size(1) != image_data.get_size(1)
                             || output_data.get_size(2) != image_data.get_size(2) || output_data.get_size(3) != image_data.get_size(3))
                    {
                        GERROR_STREAM("DicomOrientationGadget::process, inconsistent output image dimensions in array element [n s loc] = ["
                                      << n << " " << s << " " << loc << "]");
                        return GADGET_FAIL;
                    }

                    std::memcpy(&output_data(0, 0, 0, 0, n, s, loc), image_data.get_data_ptr(), image_data.get_number_of_elements() * sizeof(std::complex<float>));
                    imagearr.headers_(n, s, loc) = image_header;

                    if (has_meta && meta_ptr != NULL && imagearr.meta_.size() >= expected_headers)
                    {
                        imagearr.meta_[meta_index] = std::move(meta_container);
                    }

                    if (verbose.value())
                    {
                        GDEBUG_STREAM("DicomOrientationGadget::process(image array) stored element [n s loc] = ["
                                      << n << " " << s << " " << loc << "]");
                    }
                }
            }
        }

        if (output_initialized)
        {
            imagearr.data_ = std::move(output_data);

            if (verbose.value())
            {
                GDEBUG_STREAM("DicomOrientationGadget::process(image array) reassembled output data and is pushing downstream");
            }
        }

        if (this->next()->putq(m1) == -1)
        {
            m1->release();
            GDEBUG_STREAM("Unable to put image array on next gadgets queue");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int DicomOrientationGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
    {
        ISMRMRD::ImageHeader *img = m1->getObjectPtr();

        GadgetContainerMessage< ISMRMRD::MetaContainer >* meta = AsContainerMessage< ISMRMRD::MetaContainer >(m1->cont()->cont());

        if (this->check_correct_image_orientation_flag.value())
        {
            if (meta->getObjectPtr()->length(GADGETRON_CORRECT_IMAGE_ORIENTATION) == 0)
            {
                if (this->next()->putq(m1) == -1)
                {
                    m1->release();
                    GDEBUG_STREAM("Unable to put ismrmrd image on next gadgets queue");
                    return GADGET_FAIL;
                }

                return GADGET_OK;
            }
        }

        uint16_t data_type = img->data_type;

        if (data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            GadgetContainerMessage< hoNDArray< unsigned short > >* datamb = AsContainerMessage< hoNDArray< unsigned short > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomOrientationGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->adjust_image(m1, datamb, meta) != GADGET_OK)
            {
                GERROR_STREAM("DicomOrientationGadget::adjust_image failed for short ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            GadgetContainerMessage< hoNDArray< short > >* datamb = AsContainerMessage< hoNDArray< short > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomOrientationGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->adjust_image(m1, datamb, meta) != GADGET_OK)
            {
                GERROR_STREAM("DicomOrientationGadget::adjust_image failed for short ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_UINT)
        {
            GadgetContainerMessage< hoNDArray< unsigned int > >* datamb = AsContainerMessage< hoNDArray< unsigned int > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomOrientationGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->adjust_image(m1, datamb, meta) != GADGET_OK)
            {
                GERROR_STREAM("DicomOrientationGadget::adjust_image failed for unsigned int ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_INT)
        {
            GadgetContainerMessage< hoNDArray< int > >* datamb = AsContainerMessage< hoNDArray< int > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomOrientationGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->adjust_image(m1, datamb, meta) != GADGET_OK)
            {
                GERROR_STREAM("DicomOrientationGadget::adjust_image failed for int ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            GadgetContainerMessage< hoNDArray< float > >* datamb = AsContainerMessage< hoNDArray< float > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomOrientationGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->adjust_image(m1, datamb, meta) != GADGET_OK)
            {
                GERROR_STREAM("DicomOrientationGadget::adjust_image failed for float ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            GadgetContainerMessage< hoNDArray< double > >* datamb = AsContainerMessage< hoNDArray< double > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomOrientationGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->adjust_image(m1, datamb, meta) != GADGET_OK)
            {
                GERROR_STREAM("DicomOrientationGadget::adjust_image failed for double ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            GadgetContainerMessage< hoNDArray< std::complex<float> > >* datamb = AsContainerMessage< hoNDArray< std::complex<float> > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomOrientationGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->adjust_image(m1, datamb, meta) != GADGET_OK)
            {
                GERROR_STREAM("DicomOrientationGadget::adjust_image failed for std::complex<float> ... ");
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            GadgetContainerMessage< hoNDArray< std::complex<double> > >* datamb = AsContainerMessage< hoNDArray< std::complex<double> > >(m1->cont());
            if (!datamb)
            {
                GERROR_STREAM("DicomOrientationGadget::process, invalid image message objects");
                return GADGET_FAIL;
            }

            if (this->adjust_image(m1, datamb, meta) != GADGET_OK)
            {
                GERROR_STREAM("DicomOrientationGadget::adjust_image failed for std::complex<double> ... ");
            }
        }
        else
        {
            GERROR_STREAM("DicomOrientationGadget::process, does not supprot data type : " << data_type);
        }

        if (this->next()->putq(m1) == -1)
        {
            m1->release();
            GDEBUG_STREAM("Unable to put ismrmrd image on next gadgets queue");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(DicomOrientationGadget)

} /* namespace Gadgetron */