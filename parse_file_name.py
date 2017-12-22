
import re
"""
Basic function for parse input sample name to formatted format in order to import into pipelines.
:type Parse_Function
"""
def pfn(filename, wanted):

    filename = os.path.basename(filename)
    storged_dict = {}
    if not PROJECT_NAME:
        if '_' in filename and '-' in filename:
            storged_dict['project_name'] = filename[:min(filename.index('_'), filename.index('-'))]
        else:
            try:
                storged_dict['project_name'] = filename[:filename.index('_')]
            except:
                storged_dict['project_name'] = filename[:filename.index('-')]
    else:
        storged_dict['project_name'] = PROJECT_NAME
    try:
        storged_dict['sample_name'] = re.findall(sample_name_pattern,filename)[0]
        storged_dict['sample_ID'] = re.findall(sample_ID_pattern, filename)[0]
    except:
        print 'wrong pattern, please check it.\nOri file name: %s' % (filename)
    if 'R' in filename:
        try:
            storged_dict['pair_ID'] = re.findall(pair_ID_pattern, filename)[0]
        except:
            print 'wrong pattern, please check it.\nOri file name: %s' % (filename)
    if ((NORMAL_SIG in filename) or (TUMOR_SIG in filename)) and (NORMAL_SIG and TUMOR_SIG):
        try:
            storged_dict['mt2_for'] = re.findall(mt2_for_pattern, filename)[0]
            storged_dict['pair_name'] = re.findall(pair_name_pattern,filename)[0]
        except:
            print 'wrong pattern, please check it.\nOri file name: %s' % (filename)
    else:
        if input_pair_tuple:
            # self input mt2 pair name, for some you can't strip simple N/T out to form pair result prefix file.
            for _pair in input_pair_tuple:
                #order is important, first one is tumor sample, next is normal sample
                if _pair[0] in filename:
                    storged_dict['mt2_for'] = TUMOR_SIG
                    storged_dict['pair_name'] = storged_dict['project_name'] + '_%s' % input_pair_tuple.index(_pair)
                elif _pair[1] in filename:
                    storged_dict['mt2_for'] = NORMAL_SIG
                    storged_dict['pair_name'] = storged_dict['project_name'] + '_%s' % input_pair_tuple.index(_pair)
                else:
                    continue
    if wanted == 'all':
        return storged_dict
    else:
        return storged_dict[wanted]

from main import *


