#from sklearn import metrics
#from sklearn.cross_validation import cross_val_score,train_test_split, LeaveOneOut
#from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
import pickle as pk
###
# Load Data 
#(gene_mat,se_mat,drug_sorted,gene_sorted,se_sorted)

######### Important file to update #################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = pk.load(open('tmp/data.txt','r')) 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X 			= data['gene_mat']#[:100,:]
se_mat		= data['se_mat']#[:100,:]
drug_names	= data['drug_sorted']
se_names	= data['se_sorted']
gene_names	= data['gene_sorted']
data = []

# array lengths
n_samples,n_features = X.shape
n_models	= se_mat.shape[1]
####
#### load classifer 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#clf = joblib.load("tmp/Ascore-clf_met.out")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####
def GetSideEffects (drug_name):
  ind = drug_names.index(drug_name)
  try :
  	se_m = se_mat[ind] 
  except IndexError: 
  	return (null)
  print se_m[ind]
  return (drug_name[ind],se_m[ind])
