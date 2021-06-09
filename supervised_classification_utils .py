#!/usr/bin/env python
# coding: utf-8

# In[3]:


def classifier(data,test):
    from sklearn.preprocessing import StandardScaler
    best_columns=['s_ORF1ab', 'del_ORF1ab', 'mnc_S', 'del_S', 's_E', 's_M', 'mc_N',
       'mnc_N', 'del_N', 'mnc_NON_COD']
    
    target = data['label']
    data=data[best_columns]
    test=test[best_columns]
    
    ss=StandardScaler()
    data=ss.fit_transform(data)
    test=ss.transform(test)
    clf =RandomForestClassifier(random_state=42,max_depth=10,n_estimators=10) 
    clf.fit(data,target)
    target_test = clf.predict(test)
    return target_test
    

