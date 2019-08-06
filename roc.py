#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import scipy
import scipy.stats as stats
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import functools
import itertools
from sklearn import metrics


# print('numpy: v' + str(np.__version__))
# print('scipy: v' + str(scipy.__version__))
# print('pandas: v' + str(pd.__version__))
# print('matplotlib: v' + str(mpl.__version__))


def calculate_roc(data, tplist, tnlist=None):
    # The variable data should be a sorted series, with positive at the beginning and negative at the end.
    # The index of variable data should be the name.
    # The variable tplist and tnlist should be list or pandas.Series
    # The tplist should be gaven.
    # And if tnlist is not provided,
    # the tnlist will be set to the name from data not in the tplist as tnlist.
    fl = data.index.to_series()  # full list
    tplist = pd.Series(tplist)
    tpl = tplist[tplist.isin(fl)]  # true positve list
    tnl = []  # true negative list
    if tnlist is None:
        tnl = fl[~fl.isin(tpl)]
    else:
        tnlist = pd.Series(tnlist)
        tnl = tnlist[tnlist.isin(fl)]
    del tplist
    del tnlist

    intpl = data.index.isin(tpl)
    intnl = data.index.isin(tnl)
    fpr, tpr, thresholds = metrics.roc_curve(intpl, data, pos_label=True)

    roc_data = pd.DataFrame(
        {
            'score': thresholds,
            'tpr': tpr,
            'fpr': fpr
        }
    )
    roc_data.sort_values('fpr', ascending=True, inplace=True)
    roc_data['x_topright'] = roc_data['fpr']
    roc_data['y_topright'] = roc_data['tpr']
    roc_data['x_bottomleft'] = pd.concat(
        [pd.Series([0]), roc_data['x_topright'][0:-1]], ignore_index=True
    )
    roc_data['y_bottomleft'] = 0
    roc_data['h'] = roc_data['y_topright'] - roc_data['y_bottomleft']
    roc_data['w'] = roc_data['x_topright'] - roc_data['x_bottomleft']
    roc_data['s'] = roc_data['h'] * roc_data['w']
    roc_data.reset_index(drop=True, inplace=True)
    return roc_data


# ***

# ## Gene

def plot_roc_ess(expname, roc, esscount):
    fig, axis = plt.subplots(1, 2, figsize=(14, 8))

    # ROC
    for x in roc:
        axis[0].plot(
            1 - roc[x]['tnr'], roc[x]['tpr'],
            label='{0} AUC: {1:.4f}'.format(x, roc[x]['s'].sum())
        )

    axis[0].plot([0,1], [0,1], 'k', linestyle='dotted')
    axis[0].set_title('ROC of {}'.format(expname))
    axis[0].set_xlabel('False Positive Rate')
    axis[0].set_ylabel('True Positive Rate')
    axis[0].set_aspect(1)
    axis[0].legend()

    # Essential counts
    for x in esscount:
        axis[1].plot(
            esscount[x].index, esscount[x],
            label='{0}'.format(x)
        )

    axis[1].set_title('Known essential gene step of {}'.format(expname))
    axis[1].set_xlabel('Ranks')
    axis[1].set_ylabel('True Essential Gene Counts')
    axis[1].set_aspect(axis[1].get_xlim()[1]/axis[1].get_ylim()[1])
    axis[1].legend(loc='lower right')
    return fig


# In[4]:


dirs = dict()
dirs['data'] = 'data'
dirs['zfc'] = os.path.join(dirs['data'], 'zfc')
dirs['fig'] = 'fig'
dirs['venn'] = 'venn'

# In[5]:


ess = pd.read_table(os.path.join(dirs['data'], 'essentials.csv'), header=None)
ess.columns = ['Gene', 'ID']
noness = pd.read_table(os.path.join(dirs['data'], 'nonessentials.csv'), header=None)
noness.columns = ['Gene', 'ID']


# In[22]:


labels = [
    'HeLa_IBAR_MOI3_D15_Rm',
    'HeLa_IBAR_MOI3_D21_Rm',
    'HeLa_IBAR_MOI3_Exp_Rm',
    'HeLa_MOI3_D15_R1',
    'HeLa_MOI3_D15_R2',
    'HeLa_MOI3_D15_Rm',
    'HeLa_MOI3_D21_R1',
    'HeLa_MOI3_D21_R2',
    'HeLa_MOI3_D21_Rm',
    'HeLa_MOI3_Exp_Rm',
    'K562_MOI10_D30_Rm',
    'K562_MOI3_D30_Rm'
]


# In[23]:


label_group = {
    'HeLa_BBK': [
        'HeLa_MOI3_D15_R1',
        'HeLa_MOI3_D15_R2',
        'HeLa_MOI3_D15_Rm',
        'HeLa_MOI3_D21_R1',
        'HeLa_MOI3_D21_R2',
        'HeLa_MOI3_D21_Rm',
        'HeLa_MOI3_Exp_Rm'
    ],
    'HeLa_IBAR': [
        'HeLa_IBAR_MOI3_D15_Rm',
        'HeLa_IBAR_MOI3_D21_Rm',
        'HeLa_IBAR_MOI3_Exp_Rm',
    ],
    'HeLa_BBK_D15': [
        'HeLa_MOI3_D15_R1',
        'HeLa_MOI3_D15_R2',
        'HeLa_MOI3_D15_Rm',
    ],
    'HeLa_BBK_D21': [
        'HeLa_MOI3_D21_R1',
        'HeLa_MOI3_D21_R2',
        'HeLa_MOI3_D21_Rm',
    ],
    'HeLa_merge': [
        'HeLa_IBAR_MOI3_D15_Rm',
        'HeLa_IBAR_MOI3_D21_Rm',
        'HeLa_IBAR_MOI3_Exp_Rm',
        'HeLa_MOI3_D15_Rm',
        'HeLa_MOI3_D21_Rm',
        'HeLa_MOI3_Exp_Rm',
    ],
    'K562_BBK': [
        'K562_MOI10_D30_Rm',
        'K562_MOI3_D30_Rm'
    ]
}


# In[24]:


zfc = dict()
fc = dict()
sgrna_lfc = dict()

for x in labels:
    zfc[x] = pd.read_table(
        os.path.join(dirs['zfc'], x, '{}_gene.txt'.format(x)),
        sep='\t', header=0
    )
    fc[x] = pd.read_table(
        os.path.join(dirs['zfc'], x, '{}_sgrna.txt'.format(x)),
        sep='\t', header=0
    )
    sgrna_lfc[x] = fc[x][
        ['gene', 'guide', 'barcode', 'ctrl', 'exp']
    ].groupby(['gene', 'guide'])[['ctrl', 'exp']].sum().reset_index()
    sgrna_lfc[x]['lfc'] = np.log2(
        sgrna_lfc[x]['exp'] / sgrna_lfc[x]['ctrl']
    )


zfc_roc = dict()


for x in zfc:
    zfc_roc[x] = pd.DataFrame(
        {'score': zfc[x]['zlfc'].unique()}
    )
    zfc_roc[x] = zfc_roc[x].assign(
        tpr = zfc_roc[x]['score'].map(
            lambda a: (
                (zfc[x]['zlfc'] <= a) & (zfc[x]['gene'].isin(ess['Gene']))
            ).sum() / ess['Gene'].isin(zfc[x]['gene']).sum()
        ),
        tnr = zfc_roc[x]['score'].map(
            lambda a: (
                (zfc[x]['zlfc'] > a) & (zfc[x]['gene'].isin(noness['Gene']))
            ).sum() / noness['Gene'].isin(zfc[x]['gene']).sum()
        )
    )
    zfc_roc[x].sort_values('tnr', ascending=False, inplace=True)
    zfc_roc[x]['x_topright'] = 1 - zfc_roc[x]['tnr']
    zfc_roc[x]['y_topright'] = zfc_roc[x]['tpr']
    zfc_roc[x]['x_bottomleft'] = np.concatenate(
        [np.array([0]), zfc_roc[x]['x_topright'][0:-1]]
    )
    zfc_roc[x]['y_bottomleft'] = 0
    zfc_roc[x]['h'] = zfc_roc[x]['y_topright'] - zfc_roc[x]['y_bottomleft']
    zfc_roc[x]['w'] = zfc_roc[x]['x_topright'] - zfc_roc[x]['x_bottomleft']
    zfc_roc[x]['s'] = zfc_roc[x]['h'] * zfc_roc[x]['w']


# In[28]:


zfc_essc = dict()

for x in zfc:
    count = zfc[x].reset_index(
        drop=False
    ).sort_values(
        'zlfc', ascending=True
    )['gene'].isin(ess['Gene']).reset_index(drop=True).cumsum()
    zfc_essc[x] = count[np.arange(20, 10000, 20)]


# ## ROC

# In[29]:


for expname, plotlabels in label_group.items():

    roc_dict = dict(
        zip(
            [
                '{}'.format(x.replace('HeLa_', '').replace('K562_', ''))
                for x in plotlabels
            ],
            [
                zfc_roc[x] for x in plotlabels
            ]
        )
    )

    essc_dict = dict(
        zip(
            [
                '{}'.format(x.replace('HeLa_', '').replace('K562_', ''))
                for x in plotlabels
            ],
            [
                zfc_essc[x] for x in plotlabels
            ]
        )
    )

    fig = plot_roc_ess(expname, roc_dict, essc_dict)
    fig.savefig(
        os.path.join(dirs['fig'], '{0}_roc.pdf'.format(expname))
    )
    # fig.savefig(os.path.join(dirs['fig'], '{0}_roc.png'.format(expname)))
    # plt.show()


# ## sgRNA

# In[30]:


def calculate_cumulative_fraction(data, tplist, tnlist=None):
    # The variable data should be a sorted series, with positive at the beginning and negative at the end.
    # The index of variable data should be the name.
    # The variable tplist and tnlist should be list or pandas.Series
    # The tplist should be gaven.
    # And if tnlist is not provided,
    # the tnlist will be set to the name from data not in the tplist as tnlist.

    fl = data.index.to_series()  # full list
    tplist = pd.Series(tplist)
    tpl = tplist[tplist.isin(fl)]  # true positve list
    tnl = []  # true negative list
    if tnlist is None:
        tnl = fl[~fl.isin(tpl)]
    else:
        tnlist = pd.Series(tnlist)
        tnl = tnlist[tnlist.isin(fl)]
    del tplist
    del tnlist

    percent_rank = pd.Series(list(range(1, data.shape[0] + 1))) / data.shape[0]
    cum_sum_inlist = data.index.isin(tpl).cumsum() / data.index.isin(tpl).sum()

    curve_data = pd.DataFrame(
        {
            'x': percent_rank,
            'y': cum_sum_inlist
        }
    )
    curve_data['score'] = data.reset_index(drop=True)
    curve_data['x_topright'] = curve_data['x']
    curve_data['y_topright'] = curve_data['y']
    curve_data['x_bottomleft'] = pd.concat(
        [pd.Series([0]), curve_data['x_topright'][0:-1]], ignore_index=True
    )
    curve_data['y_bottomleft'] = 0
    curve_data['h'] = curve_data['y_topright'] - curve_data['y_bottomleft']
    curve_data['w'] = curve_data['x_topright'] - curve_data['x_bottomleft']
    curve_data['s'] = curve_data['h'] * curve_data['w']
    curve_data.reset_index(drop=True, inplace=True)
    return curve_data


# In[31]:


def plot_cumulative_fraction(expname, cp_ess, cp_noness, cp_nontarget):
    fig, axis = plt.subplots(1, 1, figsize=(14, 8))

    # ROC
    for x in cp_ess:
        axis.plot(
            cp_ess[x]['x'], cp_ess[x]['y'],
            label='{0} Essential AUC: {1:.4f}'.format(
                x, cp_ess[x]['s'].sum()
            ),
        )
        axis.plot(
            cp_noness[x]['x'], cp_noness[x]['y'],
            label='{0} Non-essential AUC: {1:.4f}'.format(
                x, cp_noness[x]['s'].sum()
            ),
            linestyle='dashed'
        )
        axis.plot(
            cp_nontarget[x]['x'], cp_nontarget[x]['y'],
            label='{0} Non-target AUC: {1:.4f}'.format(
                x, cp_nontarget[x]['s'].sum()
            ),
            linestyle='dotted'
        )

    axis.set_title('{}'.format(expname))
    axis.set_xlabel('Percent-rank of sgRNAs')
    axis.set_ylabel('Cumulative fraction')
    axis.set_aspect(1)
    axis.legend(bbox_to_anchor=(1.1, 1.0))

    return fig


# In[32]:


cp_ess = dict()
cp_noness = dict()
cp_nontarget = dict()

for x in labels:
    nontarget = sgrna_lfc[x]['gene'][
        (sgrna_lfc[x]['gene'].str.find('negCon') >= 0) |
        (sgrna_lfc[x]['gene'].str.find('negative') >= 0)
    ].unique()
    sgrna_lfc[x]['idx'] = sgrna_lfc[x]['gene'] + ':' + sgrna_lfc[x]['guide']

    data = sgrna_lfc[x][
        ['idx', 'lfc']
    ].sort_values('lfc', ascending=True).set_index('idx')['lfc']

    esslist = sgrna_lfc[x]['idx'][sgrna_lfc[x]['gene'].isin(ess['Gene'])]

    cp_ess[x] = calculate_cumulative_fraction(data, esslist)

    nonesslist = sgrna_lfc[x]['idx'][sgrna_lfc[x]['gene'].isin(noness['Gene'])]

    cp_noness[x] = calculate_cumulative_fraction(data, nonesslist)

    nontargetlist = sgrna_lfc[x]['idx'][sgrna_lfc[x]['gene'].isin(nontarget)]

    cp_nontarget[x] = calculate_cumulative_fraction(data, nontargetlist)


# ### Cumulative Fraction plot

# In[33]:


for expname, plotlabels in label_group.items():
    cpess_dict = {
        '{}'.format(x.replace('HeLa_', '').replace('K562_', '')): y
        for x,y in cp_ess.items() if x in plotlabels
    }

    cpnoness_dict = {
        '{}'.format(x.replace('HeLa_', '').replace('K562_', '')): y
        for x,y in cp_noness.items() if x in plotlabels
    }

    cpnontarget_dict = {
        '{}'.format(x.replace('HeLa_', '').replace('K562_', '')): y
        for x,y in cp_nontarget.items() if x in plotlabels
    }

    fig = plot_cumulative_fraction(
        expname, cpess_dict, cpnoness_dict, cpnontarget_dict
    )
    fig.savefig(
        os.path.join(dirs['fig'], '{0}_sgrna_cp.pdf'.format(expname))
    )
    # fig.savefig(os.path.join(dirs['fig'], '{0}_sgrna_cp.png'.format(expname)))
    # plt.show()


# ### dAUC

# In[34]:


def plot_dauc(expname, dauc_dict, auc_ess_dict, auc_noness_dict, auc_nontarget_dict):
    bw = 0.2

    plabel = list(dauc_dict.keys())
    y = list(range(len(plabel)))

    y_ess = [i + bw for i in y]
    y_noness = y
    y_nontarget = [i - bw for i in y]

    fig, axes = plt.subplots(1, 3, figsize=(8, 2 + 0.5 * len(plabel)))

    # Axes 0
    daucb = axes[0].barh(
        y, [dauc_dict[x] for x in plabel],
        height=bw * 3, color='#e41a1c'
    )
    for i, lab in enumerate(plabel):
        axes[0].annotate(
            '{:.2f}'.format(dauc_dict[lab]),
            xy=(dauc_dict[lab], y[i]),
            textcoords="offset points",
            xytext=(5, 0),
            ha='left', va='center'
        )
    axes[0].set_frame_on(False)
    axes[0].set_yticks([])
    axes[0].set_yticklabels([])
    axes[0].set_xlabel('dAUC')

    # Axes 1
    for i, pos in enumerate(y):
        axes[1].text(
            0, pos, plabel[i],
        )
    axes[1].set_ylim(axes[0].get_ylim())
    axes[1].set_frame_on(False)
    axes[1].set_yticks([])
    axes[1].set_yticklabels([])
    axes[1].set_xticks([])
    axes[1].set_xticklabels([])

    # Axes 2
    essb = axes[2].barh(
        y_ess, [auc_ess_dict[x] for x in plabel],
        height=bw, color='#377eb8',
        label='Essential'
    )
    for i, lab in enumerate(plabel):
        axes[2].annotate(
            '{:.2f}'.format(auc_ess_dict[lab]),
            xy=(auc_ess_dict[lab], y_ess[i]),
            textcoords="offset points",
            xytext=(-5, 0),
            ha='right', va='center'
        )
    nonessb = axes[2].barh(
        y_noness, [auc_noness_dict[x] for x in plabel],
        height=bw, color = '#4daf4a',
        label='Non-essential'
    )
    for i, lab in enumerate(plabel):
        axes[2].annotate(
            '{:.2f}'.format(auc_noness_dict[lab]),
            xy=(auc_noness_dict[lab], y_noness[i]),
            textcoords="offset points",
            xytext=(-5, 0),
            ha='right', va='center'
        )
    nontargetb = axes[2].barh(
        y_nontarget, [auc_nontarget_dict[x] for x in plabel],
        height=bw, color='#984ea3',
        label='Non-target'
    )
    for i, lab in enumerate(plabel):
        axes[2].annotate(
            '{:.2f}'.format(auc_nontarget_dict[lab]),
            xy=(auc_nontarget_dict[lab], y_nontarget[i]),
            textcoords="offset points",
            xytext=(-5, 0),
            ha='right', va='center'
        )
    axes[2].invert_xaxis()
    axes[2].set_yticks([])
    axes[2].set_yticklabels([])
    axes[2].set_frame_on(False)
    axes[2].set_xlabel('AUC')

    legend = fig.legend(
        [daucb, essb, nonessb, nontargetb],
        ['dAUC', 'Essential', 'Non-essential', 'Non-target'],
        ncol=4, loc='upper center'#, bbox_to_anchor=(0.5, 0.95)
    )
    legend.set_frame_on(False)
    axes[1].set_xlabel(expname)
    return fig


# In[35]:


auc_ess = dict()
auc_noness = dict()
auc_nontarget = dict()
dauc = dict()

for x in labels:
    auc_ess[x] = metrics.auc(
        cp_ess[x]['x'], cp_ess[x]['y']
    )
    auc_noness[x] = metrics.auc(
        cp_noness[x]['x'], cp_noness[x]['y']
    )
    auc_nontarget[x] = metrics.auc(
        cp_nontarget[x]['x'], cp_nontarget[x]['y']
    )
    dauc[x] = auc_ess[x] - auc_noness[x]


# In[36]:


for expname, plotlabels in label_group.items():
    dauc_dict = {
        '{}'.format(x.replace('HeLa_', '')): y
        for x, y in dauc.items() if x in plotlabels
    }
    auc_ess_dict = {
        '{}'.format(x.replace('HeLa_', '')): y
        for x, y in auc_ess.items() if x in plotlabels
    }
    auc_noness_dict = {
        '{}'.format(x.replace('HeLa_', '')): y
        for x, y in auc_noness.items() if x in plotlabels
    }
    auc_nontarget_dict = {
        '{}'.format(x.replace('HeLa_', '')): y
        for x, y in auc_nontarget.items() if x in plotlabels
    }

    fig = plot_dauc(
        expname, dauc_dict, auc_ess_dict, auc_noness_dict, auc_nontarget_dict
    )
    fig.savefig(
        os.path.join(dirs['fig'], '{0}_sgrna_cp_dauc.pdf'.format(expname))
    )
#    fig.savefig(os.path.join(dirs['fig'], '{0}_sgrna_cp_dauc.png'.format(expname)))
    # plt.show()


# In[ ]:





# ###  Venn
# using website: https://www.meta-chart.com/venn#/data
# 
# threshold: p < 0.05

# In[37]:


zfc_genes = dict()

for x in labels:
    zfc_genes[x] = zfc[x]['gene'][
        (zfc[x]['p'] < 0.05) & (zfc[x]['zlfc'] < 0)
    ]

for x in labels:
    zfc_genes[x].to_csv(
        os.path.join(dirs['venn'], '{}.txt'.format(x)),
        index=False, header=False
    )
