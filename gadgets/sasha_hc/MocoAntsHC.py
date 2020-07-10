import os
# os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = '1'
import gadgetron
import ants
import numpy as np 
import multiprocessing_on_dill as multiprocessing
import dill
import xml.etree.ElementTree as ET
from scipy.ndimage import laplace
import glob
import time
import re

from tempfile import TemporaryDirectory, mkdtemp

def apply_ants(image,transform):
    image_ant = ants.from_numpy(image)
    out_ant = ants.apply_transforms(image_ant, image_ant, transform,interpolator="bSpline")
    return out_ant.numpy()

def apply_transform(image,transform):
    output = np.zeros(image.shape, dtype=np.complex64)
    for c in range(2):
        output[c].real = apply_ants(image[c].real,transform)
        output[c].imag = apply_ants(image[c].imag,transform)
    return output


def perform_ants_registration(registration_pairs):
    outprefix = mkdtemp() + "/"
    

    args = ['-d', '2',
            '-w', '[0.01,0.99]',
            '-u', '1',
            '-t', 'BSplineSyn[0.15,4x4,0x0,1]']

    for (fixed_ants, moving_ants, metric, weight1, weight2) in registration_pairs:
        f = ants.utils.get_pointer_string(fixed_ants)
        m = ants.utils.get_pointer_string(moving_ants)
        args.extend(['-m', '%s[%s,%s,%f,%d]' % (metric, f, m, weight1, weight2)])

    args.extend(['-c', '[25x25x25,1e-3,5]'])
    args.extend(['-f', '4x2x1'])
    args.extend(['-s', '1x0.5x0'])
    # args.extend(['-v','1'])
    args.extend(['-z', '1'])
    args.extend(['--float', '1'])
    args.extend(['-o', '[%s]' % (outprefix,)])

    libfn = ants.utils.get_lib_fn("antsRegistration")
    result = libfn(args)

    afffns = glob.glob(outprefix + "*" + "[0-9]GenericAffine.mat")
    fwarpfns = glob.glob(outprefix + "*" + "[0-9]Warp.nii.gz")
    iwarpfns = glob.glob(outprefix + "*" + "[0-9]InverseWarp.nii.gz")
    if len(afffns) == 0:
        afffns = ""
    if len(fwarpfns) == 0:
        fwarpfns = ""
    if len(iwarpfns) == 0:
        iwarpfns = ""

    alltx = sorted(glob.glob(outprefix + "*" + "[0-9]*"))
    findinv = np.where(
        [re.search("[0-9]InverseWarp.nii.gz", ff) for ff in alltx]
    )[0]
    findfwd = np.where([re.search("[0-9]Warp.nii.gz", ff) for ff in alltx])[
        0
    ]
    if len(findinv) > 0:
        fwdtransforms = list(
            reversed(
                [ff for idx, ff in enumerate(alltx) if idx != findinv[0]]
            )
        )
        invtransforms = [
            ff for idx, ff in enumerate(alltx) if idx != findfwd[0]
        ]
    else:
        fwdtransforms = list(reversed(alltx))
        invtransforms = alltx

    return {
        "fwdtransforms": fwdtransforms,
        "invtransforms": invtransforms,
    }

def register_images(fixed_image, moving_image):
    registration_terms = []

    fixed_lc = ants.from_numpy(fixed_image[0])
    moving_lc = ants.from_numpy(moving_image[0])

    fixed_hc = ants.from_numpy(fixed_image[1]-fixed_image[0])
    moving_hc = ants.from_numpy(moving_image[1]-moving_image[0])

    fixed_lc_laplace = ants.iMath_laplacian(fixed_lc, sigma=1.0, normalize=True)
    moving_lc_laplace = ants.iMath_laplacian(moving_lc, sigma=1.0, normalize=True)

    fixed_hc_laplace = ants.iMath_laplacian(fixed_hc, sigma=1.0, normalize=True)
    moving_hc_laplace = ants.iMath_laplacian(moving_hc, sigma=1.0, normalize=True)

    registration_terms.append((fixed_hc_laplace, moving_hc_laplace, "Demons", 0.3, 1))
    registration_terms.append((fixed_lc_laplace, moving_lc_laplace, "Demons", 0.2, 1))

    registration_terms.append((fixed_hc, moving_hc, "CC", 0.2, 1))

    registration_terms.append((fixed_lc, moving_lc, "MI2", 0.3, 32))

    mytx = perform_ants_registration(registration_terms)

    return mytx['fwdtransforms']


def transform_data(data,transforms):
    data = np.array(data.T, order='C')
    data = np.transpose(data, [0, 2, 3, 4, 1, 5, 6])

    new_shape = data.shape

    data = np.reshape(data, (-1, data.shape[-3], data.shape[-2], data.shape[-1]))

    deformed_data = np.stack([ apply_transform(image,transform) for image,transform in zip(data,transforms)])
    deformed_data = np.reshape(deformed_data, new_shape)
    reshape_order = [0, 4, 1, 2, 3, 5, 6]
    reshape_order.reverse()
    deformed_data = np.transpose(deformed_data, reshape_order)
    return deformed_data


def user_param_to_dict(user_params):
    return {k.name: k.value_ for k in user_params}


def register_data(data):
    data = np.array(data.T, order='C')
    data = np.transpose(data, [0, 2, 3, 4, 1, 5, 6])

    new_shape = data.shape

    data = np.reshape(data, (-1, data.shape[-3], data.shape[-2], data.shape[-1]))
    adata = np.abs(data)

    deformed_data = np.array(data)

    img_base = adata[0]

    def register_single_set(n):
        transform = register_images(img_base, adata[n])
        print(transform)
        return transform



    pool = multiprocessing.Pool()
    #
    result = pool.map(register_single_set, range(data.shape[0]))
    #result = [register_single_set(n) for n in range(data.shape[0])]

    deformed_data = np.stack([ apply_transform(image,transform) for image,transform in zip(data,result)])
    deformed_data = np.reshape(deformed_data, new_shape)
    reshape_order = [0, 4, 1, 2, 3, 5, 6]
    reshape_order.reverse()
    deformed_data = np.transpose(deformed_data, reshape_order)
    return result,deformed_data


def user_param_to_dict(user_params):
    return {k.name: k.value_ for k in user_params}




def fix_metadata(metadata):
    root = ET.fromstring(metadata)
    metanode = ET.SubElement(root,'meta')
    name = ET.SubElement(metanode,'name')
    name.text = 'Skip_processing_after_recon'
    value = ET.SubElement(metanode,'value')
    value.text = 'true'

    return ET.tostring(root,encoding='ascii')

def remove_second_channel(imagearray):

    imagearray.data = imagearray.data[...,:1,:]
    imagearray.headers = imagearray.headers[:,:1,...]

    #Make sane metadata interface here. Like, with an actual N dimensional array and stuff. 
    imagearray.meta = imagearray.meta[:len(imagearray.meta)//2]

    return imagearray   


def sort_by_ending(values):
    def get_ending_number(val):
        return int(re.search('.*?([0-9]+)', val).group(1))

    return sorted(values, key=get_ending_number)


def read_user_times(header, prefix, initial_value):
    params = user_param_to_dict(header.userParameters.userParameterDouble)
    max_set = header.encoding[0].encodingLimits.set_.maximum
    max_rep = header.encoding[0].encodingLimits.repetition.maximum
    keys = sort_by_ending(filter(lambda s: s.startswith(prefix), params.keys()))
    user_times = [params[k] for k in keys]
    user_times = [initial_value] + user_times
    user_times = np.array(user_times[:max_set + 1])
    user_times = np.tile(user_times, max_rep + 1)

    return user_times

def calculate_groups(ismrmrd_header,acq_headers):
    params = user_param_to_dict(ismrmrd_header.userParameters.userParameterDouble)

    time_t2p_to_center_kspace = params['TimeT2pToCenterKspace']
    t2p_rf_duration = params['T2pRfDuration']

    max_set = ismrmrd_header.encoding[0].encodingLimits.set_.maximum
    max_rep = ismrmrd_header.encoding[0].encodingLimits.repetition.maximum

    saturation_recovery_time = read_user_times(ismrmrd_header, 'SatRecTime_', 0).ravel(order='F')

    t2_prep_duration = read_user_times(ismrmrd_header, 'T2PrepDuration', 0).ravel(order='F')

    saturation_recovery_time = correct_sat(acq_headers, saturation_recovery_time)

    group2 = np.logical_and(saturation_recovery_time == np.min(saturation_recovery_time), t2_prep_duration == 0)
    group3 = t2_prep_duration > 0

    group1 = np.logical_not(np.logical_or(group2,group3))

    return group1,group2,group3

def correct_sat(headers, sat):
    corrections = np.array([header.user_int[7] * 1e-3 for header in headers[:, 0, 0]])

    sat[corrections > 0] = corrections[corrections > 0]
    return sat

def group_registration(data,groups):
    start = time.time()

    transforms = []
    deformed = []
    for g in groups:
        transform, tmp = register_data(data[:,:,:,:,g,...])
        deformed.append(tmp)
        transforms.append(transform)


    mean_images = np.stack([np.mean(group_data,axis=4) for group_data in deformed] ,axis=4)
    mean_transform, deformed_mean = register_data(mean_images)

    deformed_data = np.zeros_like(data)
    for g,transform_list,mt in zip(groups,transforms,mean_transform):
        group_data = data[:,:,:,:,g,...]
        fixed_transforms = [t + mt for t in transform_list]
        deformed_group = transform_data(group_data,fixed_transforms)
        deformed_data[:,:,:,:,g,...] = deformed_group


    print("Time to register was ", time.time() - start)
    return deformed_data

def registration(connection):
    for image in connection:
        image_copy = image
        print("Imge shape",image.data.shape)
        groups = calculate_groups(connection.header, image.headers)
        image_copy.data = group_registration(image_copy.data,groups)
        image_copy = remove_second_channel(image_copy)
        connection.send(image_copy)
