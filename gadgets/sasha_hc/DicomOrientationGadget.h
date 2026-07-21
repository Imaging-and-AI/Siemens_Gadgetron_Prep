/** \file   DicomOrientationGadget.h
    \brief  Adjust images to standard DICOM orientation
    \author Hui Xue
*/

#ifndef DicomOrientationGadget_H
#define DicomOrientationGadget_H

#include "sashahclib_export.h"

#include <string>
#include <map>
#include <complex>
#include <iostream>
#include <mutex>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "Gadget.h"
#include "hoNDArray.h"
#include "hoNDPoint.h"
#include "ismrmrd/meta.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"

#include "ImageIOAnalyze.h"
#include "cmr_ismrmrd_util.h"

#include "mri_core_def.h"
#include "mri_core_data.h"

namespace Gadgetron
{

    /**
    Incoming images are in colume storage order
    The resulting images are in either column or row order
    Suppose a contour endo has [N 2] dimension, the first (:,0) is for the firt dimension of image and (:,1) is for the second dimension
    No matter whether image content is stored column-wise or row-wise
    */
    class EXPORTSASHAHC DicomOrientationGadget : public Gadget1Of2<ISMRMRD::ImageHeader, IsmrmrdImageArray>
    {
    public:

        typedef Gadget1Of2<ISMRMRD::ImageHeader, IsmrmrdImageArray> BaseClass;

        GADGET_DECLARE(DicomOrientationGadget);

        DicomOrientationGadget();
        virtual ~DicomOrientationGadget();

        GADGET_PROPERTY(check_correct_image_orientation_flag, bool, "Whether to check GADGETRON_CORRECT_IMAGE_ORIENTATION", false);
        GADGET_PROPERTY(output_image_in_column_storage, bool, "Whether to output imgae in colum wise order", false);
        GADGET_PROPERTY(extra_orientation_check, bool, "Whether to apply extra check for cmr views", false);
        GADGET_PROPERTY(modify_acquisition_dirs, bool, "Modify read_dir and phase_dir when adjusting orientation", false);
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

    protected:

        static void find_max_in_vec(float vector[3], float& elem, int& orientation);
        static void find_min_in_vec(float vector[3], float& elem, int& orientation);
        static void compute_norm_form_row_col(float rowvec[3], float colvec[3], float normvec[3]);
        static void extract_main_orientation(float rowvec[3], float colvec[3], float normvec[3], float& value_norm, float& value_row, float& value_column, int& orientation_main, int& orientation_row);
        // Local copy of cmr_slice_geometry::ismrmd_image_2d_to_dicom_compute_mirror_rotate.
        static void ismrmd_image_2d_to_dicom_compute_mirror_rotate(float* read_dir, float* phase_dir, float* slice_dir, float& rotate_angle, int& mirror);

        template <typename T>
        static void ismrmrd_to_world_4_corners(const ISMRMRD::ImageHeader& header, hoNDPoint<T, 3>& p1, hoNDPoint<T, 3>& p2, hoNDPoint<T, 3>& p3, hoNDPoint<T, 3>& p4)
        {
            size_t RO = header.matrix_size[0];
            size_t E1 = header.matrix_size[1];

            size_t c_ro = RO / 2;
            size_t c_e1 = E1 / 2;

            T pixel_spacing = std::max(header.field_of_view[0], header.field_of_view[1]) / std::max(header.matrix_size[0], header.matrix_size[1]);

            hoNDPoint<T, 3> pos, read_dir, phase_dir;
            pos[0] = header.position[0];
            pos[1] = header.position[1];
            pos[2] = header.position[2];
            read_dir[0] = header.read_dir[0];
            read_dir[1] = header.read_dir[1];
            read_dir[2] = header.read_dir[2];

            phase_dir[0] = header.phase_dir[0];
            phase_dir[1] = header.phase_dir[1];
            phase_dir[2] = header.phase_dir[2];

            hoNDPoint<T, 3> a, b;

            p1 = pos;
            a = read_dir;
            b = phase_dir;
            a *= (-0.5 - c_ro) * pixel_spacing;
            b *= (-0.5 - c_e1) * pixel_spacing;
            p1 += a;
            p1 += b;

            p2 = pos;
            a = read_dir;
            b = phase_dir;
            a *= (-0.5 - c_ro) * pixel_spacing;
            b *= (E1 - 1 + 0.5 - c_e1) * pixel_spacing;
            p2 += a;
            p2 += b;

            p3 = pos;
            a = read_dir;
            b = phase_dir;
            a *= (RO - 1 + 0.5 - c_ro) * pixel_spacing;
            b *= (E1 - 1 + 0.5 - c_e1) * pixel_spacing;
            p3 += a;
            p3 += b;

            p4 = pos;
            a = read_dir;
            b = phase_dir;
            a *= (RO - 1 + 0.5 - c_ro) * pixel_spacing;
            b *= (-0.5 - c_e1) * pixel_spacing;
            p4 += a;
            p4 += b;
        }

        template <typename T>
        static void world_4_corners_to_ismrmrd(ISMRMRD::ImageHeader& header, const hoNDPoint<T, 3>& p1, const hoNDPoint<T, 3>& p2, const hoNDPoint<T, 3>& p3, const hoNDPoint<T, 3>& p4)
        {
            hoNDPoint<T, 3> read_dir = p4;
            read_dir -= p1;
            read_dir.normalize();

            hoNDPoint<T, 3> phase_dir = p2;
            phase_dir -= p1;
            phase_dir.normalize();

            header.read_dir[0] = read_dir[0];
            header.read_dir[1] = read_dir[1];
            header.read_dir[2] = read_dir[2];

            header.phase_dir[0] = phase_dir[0];
            header.phase_dir[1] = phase_dir[1];
            header.phase_dir[2] = phase_dir[2];

            hoNDPoint<T, 3> pt;
            pt = p1;
            pt += p2;
            pt += p3;
            pt += p4;
            pt /= (T)(4.0);

            header.position[0] = pt[0];
            header.position[1] = pt[1];
            header.position[2] = pt[2];
        }

        template<typename T>
        static void rotate_image(float rotate_angle, ISMRMRD::ImageHeader& header, hoNDArray<T>& data, std::vector<double>& px, std::vector<double>& py)
        {
            if (std::abs(rotate_angle - 0) < 0.1f)
            {
                return;
            }

            std::vector<size_t> dim;
            data.get_dimensions(dim);

            T* pData = data.begin();

            size_t sx = data.get_size(0);
            size_t sy = data.get_size(1);

            size_t num = data.get_number_of_elements() / (dim[0] * dim[1]);
            size_t n, x, y, pt;

            size_t num_pt = px.size();
            GADGET_CHECK_THROW(num_pt == py.size());

            uint16_t matrix_size[3];
            float fov[3];

            if (std::abs(rotate_angle - 90) < 0.1f)
            {
                dim[0] = data.get_size(1);
                dim[1] = data.get_size(0);

                hoNDArray<T> res(dim);
                T* pRes = res.begin();

                for (n = 0; n < num; n++)
                {
                    for (y = 0; y < sy; y++)
                    {
                        for (x = 0; x < sx; x++)
                        {
                            pRes[n * sx * sy + sy - 1 - y + x * sy] = pData[n * sx * sy + x + y * sx];
                        }
                    }
                }

                data = res;

                hoNDPoint<double, 3> p1, p2, p3, p4;
                ismrmrd_to_world_4_corners(header, p1, p2, p3, p4);
                world_4_corners_to_ismrmrd(header, p2, p3, p4, p1);

                matrix_size[0] = header.matrix_size[1];
                matrix_size[1] = header.matrix_size[0];
                matrix_size[2] = header.matrix_size[2];

                fov[0] = header.field_of_view[1];
                fov[1] = header.field_of_view[0];
                fov[2] = header.field_of_view[2];

                if (num_pt > 0)
                {
                    for (pt = 0; pt < num_pt; pt++)
                    {
                        double rx = sy - 1 - py[pt];
                        double ry = px[pt];
                        px[pt] = rx;
                        py[pt] = ry;
                    }
                }
            }
            else if (std::abs(rotate_angle - 180) < 0.1f)
            {
                hoNDArray<T> res(dim);
                T* pRes = res.begin();

                for (n = 0; n < num; n++)
                {
                    for (y = 0; y < sy; y++)
                    {
                        for (x = 0; x < sx; x++)
                        {
                            pRes[n * sx * sy + (sx - 1 - x) + (sy - 1 - y) * sx] = pData[n * sx * sy + x + y * sx];
                        }
                    }
                }

                data = res;

                hoNDPoint<double, 3> p1, p2, p3, p4;
                ismrmrd_to_world_4_corners(header, p1, p2, p3, p4);
                world_4_corners_to_ismrmrd(header, p3, p4, p1, p2);

                matrix_size[0] = header.matrix_size[0];
                matrix_size[1] = header.matrix_size[1];
                matrix_size[2] = header.matrix_size[2];

                fov[0] = header.field_of_view[0];
                fov[1] = header.field_of_view[1];
                fov[2] = header.field_of_view[2];

                if (num_pt > 0)
                {
                    for (pt = 0; pt < num_pt; pt++)
                    {
                        double rx = sx - 1 - px[pt];
                        double ry = sy - 1 - py[pt];
                        px[pt] = rx;
                        py[pt] = ry;
                    }
                }
            }
            else if (std::abs(rotate_angle - 270) < 0.1f)
            {
                dim[0] = data.get_size(1);
                dim[1] = data.get_size(0);

                hoNDArray<T> res(dim);
                T* pRes = res.begin();

                for (n = 0; n < num; n++)
                {
                    for (y = 0; y < sy; y++)
                    {
                        for (x = 0; x < sx; x++)
                        {
                            pRes[n * sx * sy + (sx - 1 - x) * sy + y] = pData[n * sx * sy + x + y * sx];
                        }
                    }
                }

                data = res;

                hoNDPoint<double, 3> p1, p2, p3, p4;
                ismrmrd_to_world_4_corners(header, p1, p2, p3, p4);
                world_4_corners_to_ismrmrd(header, p4, p1, p2, p3);

                matrix_size[0] = header.matrix_size[1];
                matrix_size[1] = header.matrix_size[0];
                matrix_size[2] = header.matrix_size[2];

                fov[0] = header.field_of_view[1];
                fov[1] = header.field_of_view[0];
                fov[2] = header.field_of_view[2];

                if (num_pt > 0)
                {
                    for (pt = 0; pt < num_pt; pt++)
                    {
                        double rx = py[pt];
                        double ry = sx - 1 - px[pt];
                        px[pt] = rx;
                        py[pt] = ry;
                    }
                }
            }
            else
            {
                return;
            }

            compute_norm_form_row_col(header.read_dir, header.phase_dir, header.slice_dir);

            header.matrix_size[0] = matrix_size[0];
            header.matrix_size[1] = matrix_size[1];
            header.matrix_size[2] = matrix_size[2];

            header.field_of_view[0] = fov[0];
            header.field_of_view[1] = fov[1];
            header.field_of_view[2] = fov[2];
        }

        template<typename T>
        static void mirror_image(int mirror, ISMRMRD::ImageHeader& header, hoNDArray<T>& data, std::vector<double>& px, std::vector<double>& py)
        {
            if (mirror == -1) return;

            std::vector<size_t> dim;
            data.get_dimensions(dim);

            T* pData = data.begin();

            size_t sx = data.get_size(0);
            size_t sy = data.get_size(1);

            size_t num = data.get_number_of_elements() / (dim[0] * dim[1]);
            size_t n, x, y, pt;

            size_t num_pt = px.size();
            GADGET_CHECK_THROW(num_pt == py.size());

            if (mirror == 1)
            {
                hoNDArray<T> res(dim);
                T* pRes = res.begin();

                for (n = 0; n < num; n++)
                {
                    for (y = 0; y < sy; y++)
                    {
                        for (x = 0; x < sx; x++)
                        {
                            pRes[n * sx * sy + (sx - 1 - x) + y * sx] = pData[n * sx * sy + x + y * sx];
                        }
                    }
                }

                data = res;

                hoNDPoint<double, 3> p1, p2, p3, p4;
                ismrmrd_to_world_4_corners(header, p1, p2, p3, p4);
                world_4_corners_to_ismrmrd(header, p4, p3, p2, p1);

                if (num_pt > 0)
                {
                    for (pt = 0; pt < num_pt; pt++)
                    {
                        px[pt] = sx - 1 - px[pt];
                    }
                }
            }
            else if (mirror == 0)
            {
                hoNDArray<T> res(dim);
                T* pRes = res.begin();

                for (n = 0; n < num; n++)
                {
                    for (y = 0; y < sy; y++)
                    {
                        for (x = 0; x < sx; x++)
                        {
                            pRes[n * sx * sy + x + (sy - 1 - y) * sx] = pData[n * sx * sy + x + y * sx];
                        }
                    }
                }

                data = res;

                hoNDPoint<double, 3> p1, p2, p3, p4;
                ismrmrd_to_world_4_corners(header, p1, p2, p3, p4);
                world_4_corners_to_ismrmrd(header, p2, p1, p4, p3);

                if (num_pt > 0)
                {
                    for (pt = 0; pt < num_pt; pt++)
                    {
                        py[pt] = sy - 1 - py[pt];
                    }
                }
            }
            else
            {
                return;
            }

            compute_norm_form_row_col(header.read_dir, header.phase_dir, header.slice_dir);
        }

        template <typename T>
        static void image_column_row_swap(ISMRMRD::ImageHeader& mh, hoNDArray<T>& data, ISMRMRD::MetaContainer& img_attrib)
        {
            // GADGET_UNUSED(img_attrib);

            size_t sx = data.get_size(0);
            size_t sy = data.get_size(1);
            size_t num = data.get_number_of_elements() / (sx * sy);

            T* pData = data.begin();

            hoNDArray<T> buf;
            buf.create(sy, sx);
            T* pBuf = buf.begin();

            size_t x, y, n;
            for (n = 0; n < num; n++)
            {
                for (y = 0; y < sy; y++)
                {
                    for (x = 0; x < sx; x++)
                    {
                        pBuf[y + x * sy] = pData[n * sx * sy + x + y * sx];
                    }
                }

                std::memcpy(pData + n * sx * sy, pBuf, sizeof(T) * sx * sy);
            }

            std::vector<size_t> dims;
            data.get_dimensions(dims);
            dims[0] = sy;
            dims[1] = sx;
            data.reshape(dims);

            std::swap(mh.matrix_size[0], mh.matrix_size[1]);
            std::swap(mh.field_of_view[0], mh.field_of_view[1]);
        }

        virtual int process_config(ACE_Message_Block * mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);
        virtual int process(GadgetContainerMessage<IsmrmrdImageArray>* m1);

        template <typename T>
        int adjust_image_core(ISMRMRD::ImageHeader& mh, hoNDArray<T>& data, ISMRMRD::MetaContainer* img_attrib)
        {
            try
            {
                std::string filename;

                if (img_attrib != NULL)
                    create_image_name_from_header(mh, *img_attrib, filename);
                else
                    create_image_name_from_header(mh, filename);

                if (this->verbose.value())
                {
                    GDEBUG_STREAM("-->DicomOrientationGadget, incoming image : " << filename);
                }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(data, debug_folder_full_path_ + "BeforeAdjust_" + filename);
                }

                float phase_vec[3];
                phase_vec[0] = mh.phase_dir[0];
                phase_vec[1] = mh.phase_dir[1];
                phase_vec[2] = mh.phase_dir[2];

                float orig_read_dir[3];
                orig_read_dir[0] = mh.read_dir[0];
                orig_read_dir[1] = mh.read_dir[1];
                orig_read_dir[2] = mh.read_dir[2];

                float orig_phase_dir[3];
                orig_phase_dir[0] = mh.phase_dir[0];
                orig_phase_dir[1] = mh.phase_dir[1];
                orig_phase_dir[2] = mh.phase_dir[2];

                float orig_slice_dir[3];
                orig_slice_dir[0] = mh.slice_dir[0];
                orig_slice_dir[1] = mh.slice_dir[1];
                orig_slice_dir[2] = mh.slice_dir[2];

                std::vector<double> px, py;

                bool flip_image = false;

                if (!vendor_.empty())
                {
                    if (vendor_.find("Siemens") != std::string::npos
                        || vendor_.find("siemens") != std::string::npos
                        || vendor_.find("SIEMENS") != std::string::npos)
                    {
                        if (this->verbose.value())
                        {
                            GDEBUG_STREAM("Vendor - SIEMENS, flip images ... ");
                        }

                        flip_image = true;
                    }
                }

                if (flip_image)
                {
                    if (this->verbose.value())
                    {
                        GDEBUG_STREAM("Vendor - required, flip images ... ");
                    }

                    ISMRMRD::ImageHeader headerCopy(mh);
                    mirror_image(0, headerCopy, data, px, py);
                    mirror_image(1, headerCopy, data, px, py);
                }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(data, debug_folder_full_path_ + "BeforeAdjust_Flip" + filename);
                }

                float rotate_angle = 0;
                int mirror = 0;
                ismrmd_image_2d_to_dicom_compute_mirror_rotate(mh.read_dir, mh.phase_dir, mh.slice_dir, rotate_angle, mirror);

                ISMRMRD::ImageHeader mh_corr(mh);
                hoNDArray<T> data_corr(data);

                rotate_image(rotate_angle, mh_corr, data_corr, px, py);
                mirror_image(mirror, mh_corr, data_corr, px, py);

                if (extra_orientation_check.value())
                {
                    if (data_corr.get_size(0) < data_corr.get_size(1))
                    {
                        if (mh.phase_dir[0] > 0 && mh.phase_dir[1] > 0 && mh.phase_dir[2] < 0)
                        {
                            if (mh.read_dir[0] < 0 && mh.read_dir[2] < 0)
                            {
                                if (std::abs(rotate_angle - 90) < 0.1 && mirror == 1)
                                {
                                    rotate_image(90, mh_corr, data_corr, px, py);
                                    mirror_image(1, mh_corr, data_corr, px, py);
                                }
                            }
                        }
                    }
                }

                mh = mh_corr;
                data = data_corr;

                if (img_attrib != NULL)
                {
                    img_attrib->set(GADGETRON_KEEP_IMAGE_GEOMETRY, (long)1);
                }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(data, debug_folder_full_path_ + "AfterAdjust_" + filename);
                }

                if (output_image_in_column_storage.value())
                {
                    if (img_attrib != NULL)
                    {
                        image_column_row_swap(mh, data, *img_attrib);
                    }
                    else
                    {
                        ISMRMRD::MetaContainer dummy_attrib;
                        image_column_row_swap(mh, data, dummy_attrib);
                    }

                    if (!debug_folder_full_path_.empty())
                    {
                        gt_exporter_.export_array(data, debug_folder_full_path_ + "AfterAdjust_ColumnWise" + filename);
                    }
                }

                std::swap(mh.read_dir[0], mh.phase_dir[0]);
                std::swap(mh.read_dir[1], mh.phase_dir[1]);
                std::swap(mh.read_dir[2], mh.phase_dir[2]);

                compute_norm_form_row_col(mh.phase_dir, mh.read_dir, mh.slice_dir);

                if (img_attrib != NULL)
                {
                    img_attrib->set("PatientPosition", (double)mh.position[0]);
                    img_attrib->append("PatientPosition", (double)mh.position[1]);
                    img_attrib->append("PatientPosition", (double)mh.position[2]);

                    img_attrib->set("read_dir", (double)mh.read_dir[0]);
                    img_attrib->append("read_dir", (double)mh.read_dir[1]);
                    img_attrib->append("read_dir", (double)mh.read_dir[2]);

                    img_attrib->set("phase_dir", (double)mh.phase_dir[0]);
                    img_attrib->append("phase_dir", (double)mh.phase_dir[1]);
                    img_attrib->append("phase_dir", (double)mh.phase_dir[2]);

                    img_attrib->set("slice_dir", (double)mh.slice_dir[0]);
                    img_attrib->append("slice_dir", (double)mh.slice_dir[1]);
                    img_attrib->append("slice_dir", (double)mh.slice_dir[2]);

                    img_attrib->set("patient_table_position", (double)mh.patient_table_position[0]);
                    img_attrib->append("patient_table_position", (double)mh.patient_table_position[1]);
                    img_attrib->append("patient_table_position", (double)mh.patient_table_position[2]);

                    img_attrib->set("FOV", (double)mh.field_of_view[0]);
                    img_attrib->append("FOV", (double)mh.field_of_view[1]);
                    img_attrib->append("FOV", (double)mh.field_of_view[2]);

                    img_attrib->set("phase_encoding_dir", (double)phase_vec[0]);
                    img_attrib->append("phase_encoding_dir", (double)phase_vec[1]);
                    img_attrib->append("phase_encoding_dir", (double)phase_vec[2]);

                    img_attrib->set("ImageRowDir", (double)mh.read_dir[0]);
                    img_attrib->append("ImageRowDir", (double)mh.read_dir[1]);
                    img_attrib->append("ImageRowDir", (double)mh.read_dir[2]);

                    img_attrib->set("ImageColumnDir", (double)mh.phase_dir[0]);
                    img_attrib->append("ImageColumnDir", (double)mh.phase_dir[1]);
                    img_attrib->append("ImageColumnDir", (double)mh.phase_dir[2]);

                    img_attrib->set("ImageSliceDir", (double)mh.slice_dir[0]);
                    img_attrib->append("ImageSliceDir", (double)mh.slice_dir[1]);
                    img_attrib->append("ImageSliceDir", (double)mh.slice_dir[2]);

                    if (!modify_acquisition_dirs.value())
                    {
                        mh.read_dir[0] = orig_read_dir[0];
                        mh.read_dir[1] = orig_read_dir[1];
                        mh.read_dir[2] = orig_read_dir[2];

                        mh.phase_dir[0] = orig_phase_dir[0];
                        mh.phase_dir[1] = orig_phase_dir[1];
                        mh.phase_dir[2] = orig_phase_dir[2];

                        mh.slice_dir[0] = orig_slice_dir[0];
                        mh.slice_dir[1] = orig_slice_dir[1];
                        mh.slice_dir[2] = orig_slice_dir[2];
                    }
                }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array(data, debug_folder_full_path_ + "AfterAdjust_final_" + filename);
                }
            }
            catch(...)
            {
                GERROR_STREAM("Exceptions happened in DicomOrientationGadget::adjust_image_core(...) ... ");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

        template <typename T>
        int adjust_image(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< T > >* m2, GadgetContainerMessage< ISMRMRD::MetaContainer >* m3)
        {
            try
            {
                ISMRMRD::ImageHeader mh = *m1->getObjectPtr();
                hoNDArray< T > data = *m2->getObjectPtr();
                ISMRMRD::MetaContainer* img_attrib = m3 ? m3->getObjectPtr() : NULL;

                if (adjust_image_core(mh, data, img_attrib) != GADGET_OK)
                {
                    return GADGET_FAIL;
                }

                *m1->getObjectPtr() = mh;
                *m2->getObjectPtr() = data;
            }
            catch(...)
            {
                GERROR_STREAM("Exceptions happened in DicomOrientationGadget::adjust_image(...) ... ");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

    private:
        ISMRMRD::IsmrmrdHeader xml;
        std::string vendor_;
        std::string protocol_;
        std::string debug_folder_full_path_;
        Gadgetron::ImageIOAnalyze gt_exporter_;
    };

}

#endif // DicomOrientationGadget_H
