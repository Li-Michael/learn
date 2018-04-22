
# coding: utf-8

# In[1]:

get_ipython().magic(u'ls')


# In[2]:

import numpy as np
import pandas as pd
from sklearn import datasets,metrics
from sklearn.neighbors import KNeighborsClassifier


# In[3]:

# input data
total_expr = pd.read_csv('total_gene.txt',sep='\t', index_col='Probe')


# In[4]:

total_expr.iloc[:3,:3]


# In[5]:

total_t = total_expr.T


# In[6]:

total_t.iloc[:3,:3]


# In[7]:

subtotal = total_t.iloc[1:,:]
info = pd.read_csv('sampleinfo_allinfo.xls',sep='\t', index_col='sample_name')
#info.iloc[:3,8:10], subtotal.iloc[:3,:3]


# In[8]:

# exact the index that will be droped
drop_index = list(info[info['primary_site'].str.contains(r'other')].index)

# drop the index 
expr_x_src = subtotal.drop(drop_index)
expr_x = subtotal.drop(drop_index).as_matrix() # transform matrix
expr_y = info.drop(drop_index).iloc[:,8]
expr_y.iloc[:3]


# In[15]:

## PCA function
def pca(dataMat, topNfeat=999999):
    meanVals = np.mean(dataMat, axis=0)
    meanRemoved = dataMat - meanVals
    #print meanRemoved
    covMat = np.cov(meanRemoved, rowvar=0)
    eigVals, eigVects = np.linalg.eig(np.mat(covMat))
    eigValInd = np.argsort(eigVals)
    eigValInd = eigValInd[:-(topNfeat+1):-1]
    redEigVects = eigVects[:, eigValInd]
    lowDDataMat = meanRemoved * redEigVects
    reconMat = (lowDDataMat * redEigVects.T) + meanVals
    print("OK!")
    plt.plot(range(20), eigVals[:20]/sum(eigVals))
    return lowDDataMat, reconMat, eigVals, eigVects


# In[14]:

lowDMat, reconMat, eigVals, eigVects = pca(expr_x.astype(float),10)


# In[ ]:

eig


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

# split the expression data in train and test data
# a random permutation, to split the data randomly
np.random.seed(0)
indices = np.random.permutation(len(expr_x))
expr_x_train = expr_x[indices[:-10]]
expr_y_train = expr_y[indices[:-10]]
expr_x_test = expr_x[indices[-10:]]
expr_y_test = expr_y[indices[-10:]]

# creat and fit a nearest-neighbor classifier
from sklearn.neighbors import KNeighborsClassifier
knn = KNeighborsClassifier(weights='distance')
expr_x.shape, expr_y.shape
knn.fit(expr_x_train, expr_y_train)


# In[ ]:

predicted = knn.predict(expr_x_test)
predicted,expr_y_test


# In[ ]:




# In[ ]:

print("Classification report for classifier %s:\n%s\n"
      % (knn, metrics.classification_report(expr_y_test, predicted)))
print("Confusion matrix:\n%s" % metrics.confusion_matrix(expr_y_test, predicted) )


# In[ ]:




# In[ ]:



