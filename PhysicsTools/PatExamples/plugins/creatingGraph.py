from timeit import default_timer as timer
import uproot
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import math
from math import pi
import torch
import random
from torch_geometric.data import Data
import pickle
from scipy.spatial import distance
import mplhep as hep

hep.set_style(hep.style.CMS)


np.random.seed(0)


def gen_dataframe(rfilename, num_event, num_start=0):
    """
    select pfcands from original root and convert to a pandas dataframe.
    Returned is a list of dataframes, with one dataframe for one event.
    """
    print(f"reading events from {num_start} to {num_start+num_event}")
    tree = uproot.open(rfilename)["Events"]
    pfcands = tree.arrays(filter_name="PF_*", entry_start=num_start,
                          entry_stop=num_event + num_start)
    genparts = tree.arrays(filter_name="packedGenPart_*",
                           entry_start=num_start, entry_stop=num_event + num_start)
    print(tree.num_entries)

    df_pf_list = []
    df_gen_list = []

    #
    # todo: this loop can probably be removed;
    # no need to convert to dataframe for each event
    #
    for i in range(num_event):
        event = pfcands[i]
        selected_features = ['PF_eta', 'PF_phi', 'PF_pt',
                             'PF_pdgId', 'PF_charge', 'PF_puppiWeight', 'PF_puppiWeightChg', 'PF_dz',
                             'PF_fromPV'
                             ]
        pf_chosen = event[selected_features]
        df_pfcands = ak.to_dataframe(pf_chosen)
        df_pfcands = df_pfcands[abs(df_pfcands['PF_eta']) < 2.5]
        # df_pfcands['PF_pt'] = np.log(df_pfcands['PF_pt'])

        df_pf_list.append(df_pfcands)

    for i in range(num_event):
        event = genparts[i]
        selected_features = ['packedGenPart_eta', 'packedGenPart_phi',
                             'packedGenPart_pt', 'packedGenPart_pdgId', 'packedGenPart_charge']
        gen_chosen = event[selected_features]
        df_genparts = ak.to_dataframe(gen_chosen)
        # eliminate those with eta more than 2.5 and also neutrinos
        selection = (abs(df_genparts['packedGenPart_eta']) < 2.5) & (abs(df_genparts['packedGenPart_pdgId']) != 12) & (
            abs(df_genparts['packedGenPart_pdgId']) != 14) & (abs(df_genparts['packedGenPart_pdgId']) != 16)
        df_genparts = df_genparts[selection]
        df_gen_list.append(df_genparts)

    return df_pf_list, df_gen_list


def prepare_dataset(rfilename, num_event, num_start=0):
    """
    process each dataframe, prepare the ingredients for graphs (edges, node features, labels, etc).
    Returned is a list of graphs (torch.geometric data), with one graph for one event.
    """
    data_list = []
    

    df_pf_list, df_gen_list = gen_dataframe(rfilename, num_event, num_start)

    PTCUT = 0.5
    
    #defination of max edge distance
    #deltaR < 0.8/0.4 are recognized as an edge
    deltaRSetting = 0.8
    #deltaRSetting = 0.4

    for num in range(len(df_pf_list)):
        if num % 10 == 0:
            print(f"processed {num} events")
        #
        # toDO: this for loop can probably be also removed
        # and the gen_dataframe can be merged inside this function
        #
        df_pfcands = df_pf_list[num]
        # fromPV > 2 or < 1 is a really strict cut
        LV_index = np.where((df_pfcands['PF_puppiWeight'] > 0.99) & (df_pfcands['PF_charge'] != 0) & (
            df_pfcands['PF_pt'] > PTCUT) & (df_pfcands['PF_fromPV'] > 2))[0]
        PU_index = np.where((df_pfcands['PF_puppiWeight'] < 0.01) & (df_pfcands['PF_charge'] != 0) & (
            df_pfcands['PF_pt'] > PTCUT) & (df_pfcands['PF_fromPV'] < 1))[0]
        # print("LV index", LV_index)
        # print("PU index", PU_index)
        if LV_index.shape[0] < 5 or PU_index.shape[0] < 50:
            continue
        Neutral_index = np.where(df_pfcands['PF_charge'] == 0)[0]
        Charge_index = np.where(df_pfcands['PF_charge'] != 0)[0]
        # label of samples
        label = df_pfcands.loc[:, ['PF_puppiWeight']].to_numpy()
        label = torch.from_numpy(label).view(-1)
        label = label.type(torch.long)

        # node features
        node_features = df_pfcands.drop(df_pfcands.loc[:, ['PF_charge']], axis=1).drop(
            df_pfcands.loc[:, ['PF_fromPV']], axis=1).to_numpy()

        node_features = torch.from_numpy(node_features)
        node_features = node_features.type(torch.float32)

        # set the charge pdgId for one hot encoding later
        # ToDO: fix for muons and electrons
        index_pdgId = 3
        node_features[[Charge_index.tolist()], index_pdgId] = 0
        # get index for mask training and testing samples for pdgId(3) and puppiWeight(4)
        # one hot encoding for pdgId and puppiWeight
        pdgId = node_features[:, index_pdgId]
        photon_indices = (pdgId == 22)
        pdgId[photon_indices] = 1
        hadron_indices = (pdgId == 130)
        pdgId[hadron_indices] = 2
        pdgId = pdgId.type(torch.long)
        # print(pdgId)
        pdgId_one_hot = torch.nn.functional.one_hot(pdgId)
        pdgId_one_hot = pdgId_one_hot.type(torch.float32)
        assert pdgId_one_hot.shape[1] == 3, "pdgId_one_hot.shape[1] != 3"
        # print ("pdgID_one_hot", pdgId_one_hot)
        # set the neutral puppiWeight to default
        index_puppi = 4
        index_puppichg = 5
        pWeight = node_features[:, index_puppi].clone()
        pWeightchg = node_features[:, index_puppichg].clone()
        node_features[[Neutral_index.tolist()], index_puppi] = 2
        puppiWeight = node_features[:, index_puppi]
        puppiWeight = puppiWeight.type(torch.long)
        puppiWeight_one_hot = torch.nn.functional.one_hot(puppiWeight)
        puppiWeight_one_hot = puppiWeight_one_hot.type(torch.float32)
        # columnsNamesArr = df_pfcands.columns.values
        node_features = torch.cat(
            (node_features[:, 0:3], pdgId_one_hot, puppiWeight_one_hot), 1)
        #    (node_features[:, 0:3], pdgId_one_hot, node_features[:, -1:], puppiWeight_one_hot), 1)
        # i(node_features[:, 0:3], pdgId_one_hot,node_features[:,5:6], puppiWeight_one_hot), 1)
        # (node_features[:, 0:4], pdgId_one_hot, puppiWeight_one_hot), 1)

        if num == 0:
            print("pdgId dimensions: ", pdgId_one_hot.shape)
            print("puppi weights dimensions: ", puppiWeight_one_hot.shape)
            print("last dimention: ", node_features[:, -1:].shape)
            print("node_features dimension: ", node_features.shape)
            print("node_features[:, 0:3] dimention: ",
                  node_features[:, 0:3].shape)
            print("node_features dimension: ", node_features.shape)
            print("node_features[:, 6:7]",
                  node_features[:, 6:7].shape)  # dz values
            # print("columnsNamesArr", columnsNamesArr)
            # print ("pdgId_one_hot " , pdgId_one_hot)
            # print("node_features[:,-1:]",node_features[:,-1:])
            # print("puppi weights", puppiWeight_one_hot)
            # print("node features: ", node_features)

        # node_features = node_features.type(torch.float32)
        # construct edge index for graph

        phi = df_pfcands['PF_phi'].to_numpy().reshape((-1, 1))
        eta = df_pfcands['PF_eta'].to_numpy().reshape((-1, 1))

        df_gencands = df_gen_list[num]
        gen_features = df_gencands.to_numpy()
        gen_features = torch.from_numpy(gen_features)
        gen_features = gen_features.type(torch.float32)

        
        dist_phi = distance.cdist(phi, phi, 'cityblock')
        # deal with periodic feature of phi
        indices = np.where(dist_phi > pi)
        temp = np.ceil((dist_phi[indices] - pi) / (2 * pi)) * (2 * pi)
        dist_phi[indices] = dist_phi[indices] - temp
        dist_eta = distance.cdist(eta, eta, 'cityblock')

        dist = np.sqrt(dist_phi ** 2 + dist_eta ** 2)

        
        edge_source = np.where((dist < deltaRSetting) & (dist != 0))[0]
        edge_target = np.where((dist < deltaRSetting) & (dist != 0))[1]

        edge_index = np.array([edge_source, edge_target])
        edge_index = torch.from_numpy(edge_index)
        edge_index = edge_index.type(torch.long)

        graph = Data(x=node_features, edge_index=edge_index, y=label)
        graph.LV_index = LV_index
        graph.PU_index = PU_index
        graph.Neutral_index = Neutral_index
        graph.Charge_index = Charge_index
        graph.num_classes = 2
        graph.GenPart_nump = gen_features
        graph.pWeight = pWeight
        graph.pWeightchg = pWeightchg
        data_list.append(graph)

    return data_list


def main():
    start = timer()

    iname = "Wjets_output_10.root"
    num_events_train = 4000
    oname = "../data_pickle/dataset_graph_puppi_Wdebugging" + str(num_events_train)
    dataset_train = prepare_dataset(iname, 4000, 24000)
    #ptlist = np.array(dataset_train.x[:, 2])
    chglvpt = []
    chgpupt = []
    neupt = []
    chglveta = []
    chgpueta = []
    neueta = []
    chglvphi = []
    chgpuphi = []
    neuphi = []
    for j in range(len(dataset_train)):
        ptlist = np.array(dataset_train[j].x[:, 2])
        etalist = np.array(dataset_train[j].x[:, 0])
        philist = np.array(dataset_train[j].x[:, 1])
        charge_index = dataset_train[j].Charge_index
        neutral_index = dataset_train[j].Neutral_index
        lv_index = dataset_train[j].LV_index
        pu_index = dataset_train[j].PU_index
        chglv_index = list(set(lv_index) & set(charge_index))
        chgpu_index = list(set(pu_index) & set(charge_index))
        for p in chglv_index:
            chglvpt.append(ptlist[p])
            chglveta.append(etalist[p])
            chglvphi.append(philist[p])
        for p in chgpu_index:
            chgpupt.append(ptlist[p])
            chgpueta.append(etalist[p])
            chgpuphi.append(philist[p])
        for p in neutral_index:
            neupt.append(ptlist[p])
            neueta.append(etalist[p])
            neuphi.append(philist[p])

    linewidth = 1.5
    fontsize = 18

   #  %matplotlib inline
    plt.style.use(hep.style.ROOT)
    fig = plt.figure(figsize=(10, 8))
    pthist = np.array(chglvpt)
    n1, bins,patches = plt.hist(pthist, bins=40, range=(0, 40), histtype='step', color='blue', linewidth=linewidth,
             density=False, label=r'Charged LV particle counts:'+str(len(pthist)))
    pthist = np.array(chgpupt)
    n2, bins,patches = plt.hist(pthist, bins=40, range=(0, 40), histtype='step', color='red', linewidth=linewidth, 
             density=False, label=r'Charged PU particle counts:'+str(len(pthist)))
    pthist = np.array(neupt)
    n3, bins,patches= plt.hist(pthist, bins=40, range=(0, 40), histtype='step', color='green', linewidth=linewidth, 
             density=False, label=r'Neutral particle counts:'+str(len(pthist)))
    bincentres = bins[:-1]+0.5*1
    fig = plt.figure(figsize=(10, 8))
    pthist = np.array(chglvpt)
    plt.hist(bincentres, bins=bins, weights=n1/num_events_train, label=r'Charged LV particle counts:'+str(len(pthist)),histtype="step", color='blue', linewidth=linewidth)
    pthist = np.array(chgpupt)
    plt.hist(bincentres, bins=bins, weights=n2/num_events_train, label=r'Charged PU particle counts:'+str(len(pthist)),histtype="step", color='red', linewidth=linewidth)
    pthist = np.array(neupt)
    plt.hist(bincentres, bins=bins, weights=n3/num_events_train, label=r'Neutral particle counts:'+str(len(pthist)),histtype="step", color='green', linewidth=linewidth)
    plt.hist
    plt.xlabel(r"pT [GeV]")
    plt.ylabel('counts/EvtNum')
    plt.rc('legend', fontsize=fontsize)
    plt.legend()
    plt.savefig("Training_pt.pdf")
    plt.show()

    fig = plt.figure(figsize=(10, 8))
    etahist = np.array(chglveta)
    n1, bins,patches = plt.hist(etahist, bins=40, range=(-3, 3), histtype='step', color='blue', linewidth=linewidth,
             density=False, label=r'Charged LV particle counts:'+str(len(etahist)))
    etahist = np.array(chgpueta)
    n2, bins,patches = plt.hist(etahist, bins=40, range=(-3, 3), histtype='step', color='red', linewidth=linewidth, 
             density=False, label=r'Charged PU particle counts:'+str(len(etahist)))
    etahist = np.array(neueta)
    n3, bins,patches= plt.hist(etahist, bins=40, range=(-3, 3), histtype='step', color='green', linewidth=linewidth, 
             density=False, label=r'Neutral particle counts:'+str(len(etahist)))
    bincentres = bins[:-1]+0.5*6/40
    fig = plt.figure(figsize=(10, 8))
    etahist = np.array(chglveta)
    plt.hist(bincentres, bins=bins, weights=n1/num_events_train, label=r'Charged LV particle counts:'+str(len(etahist)),histtype="step", color='blue', linewidth=linewidth)
    etahist = np.array(chgpueta)
    plt.hist(bincentres, bins=bins, weights=n2/num_events_train, label=r'Charged PU particle counts:'+str(len(etahist)),histtype="step", color='red', linewidth=linewidth)
    etahist = np.array(neueta)
    plt.hist(bincentres, bins=bins, weights=n3/num_events_train, label=r'Neutral particle counts:'+str(len(etahist)),histtype="step", color='green', linewidth=linewidth)
    plt.hist
    plt.xlabel(r"eta")
    plt.ylabel('counts/EvtNum')
    plt.rc('legend', fontsize=fontsize)
    plt.legend()
    plt.savefig("Training_eta.pdf")
    plt.show()

    fig = plt.figure(figsize=(10, 8))
    phihist = np.array(chglvphi)
    n1, bins,patches = plt.hist(phihist, bins=40, range=(-3, 3), histtype='step', color='blue', linewidth=linewidth,
             density=False, label=r'Charged LV particle counts:'+str(len(phihist)))
    phihist = np.array(chgpuphi)
    n2, bins,patches = plt.hist(phihist, bins=40, range=(-3, 3), histtype='step', color='red', linewidth=linewidth, 
             density=False, label=r'Charged PU particle counts:'+str(len(phihist)))
    phihist = np.array(neuphi)
    n3, bins,patches= plt.hist(phihist, bins=40, range=(-3, 3), histtype='step', color='green', linewidth=linewidth, 
             density=False, label=r'Neutral particle counts:'+str(len(phihist)))
    bincentres = bins[:-1]+0.5*6/40
    fig = plt.figure(figsize=(10, 8))
    phihist = np.array(chglvphi)
    plt.hist(bincentres, bins=bins, weights=n1/num_events_train, label=r'Charged LV particle counts:'+str(len(phihist)),histtype="step", color='blue', linewidth=linewidth)
    phihist = np.array(chgpuphi)
    plt.hist(bincentres, bins=bins, weights=n2/num_events_train, label=r'Charged PU particle counts:'+str(len(phihist)),histtype="step", color='red', linewidth=linewidth)
    phihist = np.array(neuphi)
    plt.hist(bincentres, bins=bins, weights=n3/num_events_train, label=r'Neutral particle counts:'+str(len(phihist)),histtype="step", color='green', linewidth=linewidth)
    plt.hist
    plt.xlabel(r"phi")
    plt.ylabel('counts/EvtNum')
    plt.rc('legend', fontsize=fontsize)
    plt.legend()
    plt.savefig("Training_phi.pdf")
    plt.show()

    
    # save outputs in pickle format
    with open(oname, "wb") as fp:
        pickle.dump(dataset_train, fp)

    """num_events_test = 4000
    oname = "../data_pickle/dataset_graph_puppi_test_WjetsDR8" + str(num_events_test)
    dataset_test = prepare_dataset(iname, num_events_test, num_events_train)
    with open(oname, "wb") as fp:
        pickle.dump(dataset_test, fp)

    num_events_valid = 4000
    oname = "../data_pickle/dataset_graph_puppi_val_WjetsDR8" + str(num_events_valid)
    dataset_valid = prepare_dataset(
        iname, num_events_valid, num_events_train + num_events_test)
    with open(oname, "wb") as fp:
        pickle.dump(dataset_valid, fp)

    end = timer()
    program_time = end - start
    print("generating graph time " + str(program_time))"""


if __name__ == '__main__':
    main()
