#!/usr/bin/env python
# coding: utf-8

# In[1]:


def classifier(data,test):
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestClassifier
    import seaborn as sns
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    best_columns=['s_ORF1ab', 'del_ORF1ab', 'mnc_S', 'del_S', 's_E', 's_M', 'mc_N',
       'mnc_N', 'del_N', 'mnc_NON_COD']#Best columns found with RFE,see also Classifiers Folder
    control=data.copy()
    data2=pd.concat([control,test])
    data2.drop(columns='label',inplace=True)
    target = data['label']
    data=data[best_columns]
    test=test[best_columns]
    
    ss=StandardScaler()#standardize data
    data=ss.fit_transform(data)
    test=ss.transform(test)
    clf =RandomForestClassifier(random_state=42,max_depth=10,n_estimators=10) #fit the random forest
    clf.fit(data,target)
    target_test = clf.predict(test)#make predictions
    
    map = {0 : "Original", 1 : "Californian", 2 : "Brazilian", 3 : "English", 4 : "Nigerian", 5 : "South African"}
    pred=[]
    for line in target_test:
        pred.append(map[line])
    pred_df=pd.DataFrame()#create the output dataframe with the correspondi label for each line of the given dataset
    pred_df['Variant']=pred 
    plot=sns.countplot(x='Variant',data=pred_df,palette='Paired')
    fig = plot.get_figure()
    fig.savefig('./Output/classification_report_plot.png')
    
    data2['Predicted']= np.concatenate((control['label'],pred))
    label=data2.pop('Predicted')
    label.astype(np.int32)
    
    gene_data=data2.groupby((np.arange(len(data2.columns)) // 7) + 1, axis=1).sum()
    
    col_names=["ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10",'NON_COD']
    gene_data.columns=col_names

    gene_data['label']=label
    g = []
    for i in range(6):
      gene = gene_data[gene_data['label']==i].mean()
      g.append(gene)
    g = pd.DataFrame(g)

    lis = []
    map = {0 : "Original", 1 : "Californian", 2 : "Brazilian", 3 : "English", 4 : "Nigerian", 5 : "South African"}
    
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
    
    return pred_df
    

