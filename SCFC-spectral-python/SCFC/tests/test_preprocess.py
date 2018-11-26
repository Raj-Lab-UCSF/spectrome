#all the testing functions, excuse me for the mess!

# labelfile = 'desikan_atlas_68.csv'
# label_path = pth.get_sibling_path('dictionaries')
# label_filename = os.path.join(label_path, labelfile)
# regions, coords = get_desikan(label_filename)
# print(regions)
# print(coords)

# labelfile = 'mean80_fibercount.csv'
# outlist = 'HCP_list.h5'
# label_path = pth.get_sibling_path('data')
# label_filename = os.path.join(label_path, labelfile)
# out_filename = os.path.join(label_path, outlist)
# regions = get_HCP_order(label_filename, save = True, fileout = out_filename)
# print(regions)

# outlist = 'HCP_list.h5'
# path = pth.get_sibling_path('dictionaries')
# datapath = pth.get_sibling_path('data')
# confile = os.path.join(datapath, 'mean80_fibercount.csv')
# distfile = os.path.join(datapath, 'mean80_fiberlength.csv')
# filename = os.path.join(path, outlist)
# reorder_connectome(confile, distfile)



# '''example of the use of this -- preprocessed Chang's data accordingly'''
# datapath = '/Users/Megan/RajLab/MEG-chang'
# directories = pth.walk_tree(datapath)
# coord_filename = 'DK_coords_meg.mat'
# data_filename = 'DK_timecourse_20.mat'
# out_coords = 'DK_coords_meg.h5'
# out_data = 'DK_timecourse_20.h5'
#
# labelfile = 'OrderingAlphabetical_68ROIs.txt'
# label_path = pth.get_sibling_path('dictionaries')
# label_filename = os.path.join(label_path, labelfile)
#
# for dir in directories:
#     abspath = os.path.join(datapath,dir)
#     coord_path = os.path.join(abspath, coord_filename)
#     data_path = os.path.join(abspath, data_filename)
#
#     data_dict = add_key_data(label_filename, data_path)
#     pth.save_hdf5(os.path.join(abspath, out_data), data_dict)
#
#     coord_dict = add_key_coords(label_filename, coord_path)
#     pth.save_hdf5(os.path.join(abspath, out_coords), coord_dict)
