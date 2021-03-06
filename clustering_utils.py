#!/usr/bin/env python
# coding: utf-8

# In[1]:

#Function to compute the number of null values for each columns of the given dataset
def nullscan(df_check):
    df_nulls = (df_check == 0) 
    nulls_per_col = df_nulls.sum(axis=0)
    nulls_per_col /= len(df_check.index)
    for line,col_name in zip(nulls_per_col,df_check.columns):
        if line>0.99:#If the number of null values is >99% delete the columns
            df_check.drop(col_name,axis=1,inplace=True)
    return df_check


#Function to perform clustering on the dataset
def clustering(control,test):
    import pandas as pd
    import numpy as np
    import seaborn as sns
    from sklearn.neighbors import NearestNeighbors
    from sklearn.cluster import DBSCAN
    from sklearn import metrics
    from statistics import mode
    import matplotlib.pyplot as plt
    #Concat the control dataset with the one given by the user(or in case the default ones)
    total=pd.concat([control,test])
    data2=total.drop('label',axis=1,inplace=False)#Drop the true labels
    total=nullscan(total)#Check null values
    data=total.drop('label',axis=1,inplace=False)#Drop the true labels
    
    label=total['label']
    min_points=20#Set min points values for DBSCAN
    neighbors = NearestNeighbors(n_neighbors=2*len(total.columns),metric='euclidean')#Calculate the distances from the k nearest neighbors for each data 
    neighbors_fit = neighbors.fit(data) 
    distances, indices = neighbors_fit.kneighbors(data) #calculate the distances from the k-th neighbohrs
    distances = np.sort(distances, axis=0)
    distances = distances[:,1]
    
    
    
    
    dist=pd.Series(distances)

    der1 = dist.diff()#calculate the consecutive difference of each distance (a sort of first derivative)

    th=30 #number of zero values set as threshold 
    counter = 0
    #Procedure to estimate the best eps number(well described in report)
    for i,k in enumerate(der1[:len(data)-1]):
        if k > -1:
            if k == 0:#if k=0 count the number of consecutive zero(step length)

                counter = counter + 1
            else:
                if counter < th: # if the number of consecutive zero is less the the th
                    best_eps=distances[i-1]#update the best epsilon found
                    break
                else:
                    counter = 0 #go to the next sequence of zero
        
    af = DBSCAN(eps=best_eps,min_samples=min_points,metric='euclidean').fit(data)#perform DBSCAN
    
    labels = af.labels_
    n_clusters_ = len(np.unique(labels))
    from sklearn import metrics
    f=open('./Output/clustering_performance.txt',"w+")
    f.write('Estimated number of clusters: %d \n' % n_clusters_)
    f.write("Homogeneity: %0.3f \n" % metrics.homogeneity_score(label, labels))
    f.write("Completeness: %0.3f \n" % metrics.completeness_score(label, labels))
    f.write("V-measure: %0.3f \n" % metrics.v_measure_score(label, labels))
    f.write("Adjusted Rand Index: %0.3f \n"
          % metrics.adjusted_rand_score(label, labels))
    f.write("Adjusted Mutual Information: %0.3f \n"
          % metrics.adjusted_mutual_info_score(label, labels))
    f.write("Silhouette Coefficient %0.3f \n"
          % metrics.silhouette_score(data, labels, metric='euclidean'))
    f.close()
    
    n_variants=len(pd.unique(control.label)) #number of known variants 
    n_clusters = len(pd.unique(labels))# number of cluster founded 
    cl_dict={}
    
    possible_clusters=np.sort(pd.Series(labels).unique())
    cluster_out = []
    for l in possible_clusters:
        cl = label[labels==l]
        cluster_out.append(cl)
    
    for i in range(1,n_clusters): #for each cluster found with DBSCAN algo, assigne the label corresponding to the mode of the true label in the test set 
        if mode(cluster_out[i])==-1:#if the mode correspond to the unlabled data (-1) a new variant is found
            cl_dict[i-1]=n_variants #update the total num variants
            n_variants=n_variants+1
        else:
            cl_dict[i-1]=mode(cluster_out[i]) #else assigne the class mode to the cluster
    cl_dict[-1]=-1
    
    prediction = []
    for i in range(0,len(test)):
        prediction.append((cl_dict[labels[len(control)+i]]))  #Store the prediction in an array 
    prediction_pair=pd.DataFrame({'true_label':label,'pred':labels})#compare true label with the one founded by clustering algo
    
   
    
    # Create crosstab: ct
    ct = pd.crosstab(prediction_pair['true_label'], prediction_pair['pred'])

    #compute the confusion matrix
    confusion_matrix=np.array(ct)

    m=[]
    #visualize result 
    map = {0 : "Original", 1 : "Epsilon", 2 : "Gamma", 3 : "Alpha", 4 : "Eta", 5 : "Beta",-1:"Test samples"}
    n_new_variants = n_variants-len(pd.unique(control.label))
    map_cluster = {-1:"Outliers",0 : "Original", 1 : "Epsilon", 2 : "Gamma", 3 : "Alpha", 4 : "Eta", 5 : "Beta"}
    for i in range(n_new_variants):
        map_cluster[len(pd.unique(control.label))+i] = 'New variant ' + str(i+1)
    for variant in range(-1,6):
      for cluster in range(-1,n_clusters-1):
        m.append([confusion_matrix[variant+1,cluster+1],map[variant],map_cluster[cluster]])
    df=pd.DataFrame(m,columns=["Value", "Variant", "Cluster"])

    c = sns.color_palette('Paired')
    plot = sns.catplot(x="Cluster", y="Value", hue="Variant", kind="bar", data=df,height=7, aspect=1.4,palette="Paired")
    plt.ylabel('Number of samples', fontsize=13)
    plt.xlabel('Cluster', fontsize=13)
   
    fig = plot.fig
    fig.savefig('./Output/clustering_report_plot.png', bbox_inches='tight')
    
    
    
    
    data2['Predicted']= np.concatenate((control['label'],prediction))
    label=data2.pop('Predicted')
    label.astype(np.int32)
    
    gene_data=data2.groupby((np.arange(len(data2.columns)) // 7) + 1, axis=1).sum()
    
    col_names=["ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10",'NON_COD']
    gene_data.columns=col_names

    gene_data['label']=label
    g = []
    for i in range(n_variants):
      gene = gene_data[gene_data['label']==i].mean()
      g.append(gene)
    g = pd.DataFrame(g)

    lis = []
    map = {0 : "Original", 1 : "Epsilon", 2 : "Gamma", 3 : "Alpha", 4 : "Eta", 5 : "Beta"}
    for i in range(n_new_variants):
        map[len(pd.unique(control.label))+i] = 'New variant ' + str(i+1)
    g.drop(columns='label',inplace=True)
    for col in g:
      for n, el in enumerate(g[col]):
        row = [el,col,map[n]]
        lis.append(row)
    g2 = pd.DataFrame(lis, columns=["Number","Gene","Variant"])
    plot2 = sns.catplot(x="Gene", y="Number", hue="Variant", kind="bar", data=g2,height=7, aspect=1.4,palette="Paired")
    plt.ylabel('Avg number of mutations', fontsize=13)
    plt.xlabel('Gene', fontsize=13)
    plt.xticks(rotation = 45)
    fig2 = plot2.fig
    fig2.savefig('./Output/mutation_report_plot.png', bbox_inches='tight')
    
    return [map_cluster[p] for p in prediction], n_new_variants

