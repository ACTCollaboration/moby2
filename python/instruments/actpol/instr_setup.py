from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
pointing_fields = [
    {'name': 'az',  # some 2007 tods don't have Enc_Az_Deg_Astro
     'field': 'Enc_Az_Deg',
     'dtype': np.float64,
     'scale': np.pi / 180, 
     'shift': np.pi,
     },
    {'name': 'alt',
     'field': 'Enc_El_Deg',
     'dtype': np.float64,
     'scale': np.pi / 180, 
     },
    {'name': 'ctime',
     'field': 'C_Time',
     'dtype': np.float64,
     },
    {'name': 'enc_flags',
     'field': 'enc_flags',
     'dtype': np.uint16,
     },
    {'name': 'data_mask',
     'field': 'tes_flags',
     'dtype': np.int32,
     },
]

AR_NAMES = {
            'mbac145': 'AR1',
            'mbac215': 'AR2',
            'mbac280': 'AR3',
            'ar1': 'actpol1',
            'ar2': 'actpol2',
}

ACT_TIME_CONSTANTS = {
    'season2007': {
        'mbac145': ('2007/tcnst_20071027_row.dat', True),
        },
    'season2008': {
        'mbac145': ('2008/ar1_tau_090423.txt', False),
        'mbac215': ('2008/ar2_tau_090423.txt', False),
        'mbac280': ('2008/timeconst_ar3_2008_090802.txt', True),
        },
    'season2009': {
        'mbac145': ('2009/template_ar1_2009_090811.txt.tau', False),
        'mbac215': ('2009/template_ar2_2009_090811.txt.tau', False),

        'mbac280': ('2009/timeconst_ar3_2009_090825.txt', True),
        },
    'season2010': {
        'mbac145': ('2010/template_ar1_2010_101017.txt.tau', False),
        'mbac215': ('2010/template_ar2_2010_101017.txt.tau', False),
        'mbac280': ('2010/template_ar3_2010_100920.txt.tau', False),
        },
}

ACT_FLAT_FIELDS = {
    'season2007': {
        'mbac145': '2007/ff_mbac145_2007_v1.dict',
        },
    'season2008': {
        'mbac145': '2008/ff_mbac145_2008_v1.dict',
        'mbac215': '2008/ff_mbac215_2008_v1.dict',
        'mbac280': '2008/ff_mbac280_2008_v0.dict',
        },
    'season2009': {
        'mbac145': '2009/ff_mbac145_2009_v1.dict',
        'mbac215': '2009/ff_mbac215_2009_v1.dict',
        'mbac280': '2008/ff_mbac280_2008_v0.dict',
        },
    'season2010': {
        'mbac145': '2010/ff_mbac145_2010_v0.dict',
        'mbac215': '2010/ff_mbac215_2010_v0.dict',
        'mbac280': '2008/ff_mbac280_2008_v0.dict',
        },
    '2013': {
        'ar1': '2013/ff_ar1_2013_null.dict',
        },
    '2014': {
        'ar1': '2014/ff_ar1_2014_null.dict',
        'ar2': '2014/ff_ar2_2014_null.dict',
        },
}

ACT_SYNC_PARAMS = {
    'mbac145': 'syncDefault145.par',
    'mbac215': 'syncDefault215.par',
    'mbac180': 'syncDefault280.par',
}


