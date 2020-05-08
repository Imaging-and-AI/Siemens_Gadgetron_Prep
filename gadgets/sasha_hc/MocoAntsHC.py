import os
os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = '1'
import gadgetron
import pickle
import ants
import numpy as np 
import multiprocessing_on_dill as multiprocessing
import dill
import xml.etree.ElementTree as ET

def register_images(fixed_image,moving_image):
   
    fixed_ant = ants.from_numpy(fixed_image[0])
    moving_ant = ants.from_numpy(moving_image[0])

    fixed_ant2 = ants.from_numpy(fixed_image[1])
    moving_ant2 = ants.from_numpy(moving_image[1])

    mytx = ants.registration(fixed_ant,moving_ant,type_of_transform='SyNOnly',syn_sampling=32,syn_metric='MI',reg_iterations=[400,200,100,40,40],flow_sigma=0,total_sigma=2.0,grad_step=0.8,
        multivariate_extras=[("MI",fixed_ant2,moving_ant2,1.0,32)])
    
    def apply_ants(image):
        image_ant = ants.from_numpy(image)
        out_ant = ants.apply_transforms(fixed_ant,image_ant,mytx['fwdtransforms'])
        return out_ant.numpy()

    return apply_ants

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

def registration(connection):
    for image in connection:
        image_copy = image
        image_copy.data = register_data(image.data)
        image_copy = remove_second_channel(image_copy)
        connection.send(image_copy)

def register_data(data):
    data = np.array(data.T,order='C')
    data = np.transpose(data, [0,2,3,4,1,5,6])

    new_shape = data.shape

    data = np.reshape(data,(-1, data.shape[-3],data.shape[-2],data.shape[-1]))
    adata = np.abs(data)

    deformed_data = np.array(data)

    img_base = adata[0]

    def register_single_set(n):
        transform = register_images(img_base,adata[n])
        output = np.zeros(img_base.shape,dtype=np.complex64)
        for c in range(2):
            output[c].real = transform(data[n,c].real)
            output[c].imag = transform(data[n,c].imag)

        return output

    pool = multiprocessing.Pool()

    result = pool.map(register_single_set,range(1,data.shape[0]))
    result = [img_base] + result

    deformed_data = np.stack(result)
    deformed_data = np.reshape(deformed_data,new_shape)
    reshape_order = [0,4,1,2,3,5,6]
    reshape_order.reverse()
    deformed_data = np.transpose(deformed_data,reshape_order)
    return deformed_data

 