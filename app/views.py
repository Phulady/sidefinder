from flask import Flask,render_template,request
import sys
#from app import app
import numpy as np
import pymysql as mdb
import matplotlib.pyplot as plt
import time
import pickle as pk
from matplotlib_venn import venn2, venn2_circles
import json
#import GetPrediction as pre
#from sklearn.externals import joblib

app  = Flask(__name__)

# LOAD predictions
#pred_file = '../../008_BuildReg5FoldModel/SavedPredict/drugs_predSE_LogReg.out'
pred_file = '../../008_BuildReg5FoldModel/SavedPredict/drugs_predSE_RanFor.out'

######### for later use: need to added sco
#score_file = '../../008_BuildReg5FoldModel/SavedPredict/Ascore_SE1000_LogReg5fold.out' #
score_file = '../../008_BuildReg5FoldModel/SavedPredict/Ascore_SE1000_RanFor5fold.out'
###################
#TOP_predicted_se = []
#TOP_Y_probability = []


#db = mdb.connect(user="root", host="localhost", passwd='toor', db="SideEffectD1", charset='utf8')
##########################
def GetDrugList ():
  db = mdb.connect(user="root", host="localhost", passwd='toor', db="SideEffectD1", charset='utf8')
  with db:
	  	cur = db.cursor()
	  	cur.execute("select * from Drugs where se_count != 0") #'%s'" %drug_id)
	  	query_results = cur.fetchall()
  drug_list = []
  for result in query_results:
	  	#print result
	  	drug_list.append(dict(drug_id=result[0],drug_name=result[1],gene_count=result[2],se_count=[3]))
  db.close()
  return (drug_list)

#########################
def GetDrugName (drug_id):
  db = mdb.connect(user="root", host="localhost", passwd='toor', db="SideEffectD1", charset='utf8')
  with db:
      cur = db.cursor()
      cur.execute("select * from Drugs where drug_id = '%s'" %drug_id)
      query_results = cur.fetchall()
  drug_list = []
  for result in query_results:
      #print result
      drug_list.append(dict(drug_id=result[0],drug_name=result[1],gene_count=result[2],se_count=[3]))
  db.close()
  return (drug_list)
##################
def FromDrugListGetNameList (drug_list,drug_ids):
  drug_name_dict={}
  for rec in drug_list:
    drug_name_dict[rec['drug_id']] = rec['drug_name'];
  drug_names = [drug_name_dict[id] for id in drug_ids] 
  return(drug_names)
#########################
def GetDrugMatches (drug_id,limit):
  db = mdb.connect(user="root", host="localhost", passwd='toor', db="SideEffectD1", charset='utf8')
  with db:
      cur = db.cursor()
      # select only drugs with se available 
      cur.execute("""select drug_id2, match_count from Matches left join Drugs as d on (d.drug_id = Matches.drug_id2) where d.se_count != 0 and drug_id1 = %s ORDER BY match_count DESC LIMIT %s""", (drug_id,limit))
      query_results = cur.fetchall()
  match_list = []
  match_drug_ids = []
  for result in query_results:
      #print result
      match_list.append(result[1])#dict(drug_id2=result[0],match_count=result[1]))
      match_drug_ids.append(result[0])
  db.close()
  return (match_list,match_drug_ids)

#########################
def GetSEforMatched (match_drug_ids):
  db = mdb.connect(user="root", host="localhost", passwd='toor', db="SideEffectD1", charset='utf8')
  with db:
      cur = db.cursor()
      drugs_str = ",".join([str(i) for i in match_drug_ids])
      #print drugs_str
      cur.execute("""select drug_id, group_concat(se SEPARATOR '|') from SideEffectsNor  where drug_id in ("""+drugs_str+""") group by drug_id""")
      query_results = cur.fetchall()
  se_for_matched = []
  #print query_results
  for result in query_results:
      #print result
      se_for_matched.append(dict(drug_id=result[0],se=result[1]))#,se_count=result[2]))
  db.close()
  return (se_for_matched)

#########################
def GetSEforForTheDrug (drug_id):
  db = mdb.connect(user="root", host="localhost", passwd='toor', db="SideEffectD1", charset='utf8')
  with db:
      cur = db.cursor()
      cur.execute("select se from SideEffectsNor  where drug_id = '%s'" %drug_id)
      query_results = cur.fetchall()
  known_se_names = []
  #print query_results
  for result in query_results:
      #print result
      known_se_names.append(result[0])#,se_count=result[2]))
  db.close()
  return (known_se_names)
############################
def GetSE_mat(match_drug_ids,se_dict,se_set_arr,no_se,no_dg):
  se_mat  = np.zeros((no_dg,no_se),dtype=np.int16)
  #print se_mat,
  # Mask and put 1 for present side effects
  i=0
  for drug_id in match_drug_ids: #range(no_dg):
    mask = np.in1d(se_set_arr,np.array(se_dict[drug_id]))
    se_mat[i][mask] = 1
    i=i+1
  # sort indexes
  sum_se = se_mat.sum(axis=0)
  se_ind_sort = sorted(range(len(sum_se)),key=lambda x:sum_se[x],reverse=True) # get indexes of order
  print "Shape of se_mat:", se_mat.shape
  print "sum_se, se_ind_sort:",sum_se,se_ind_sort
  se_mat_sorted = se_mat[:,se_ind_sort]
  return (se_mat_sorted,se_ind_sort)

#############################
def PlotFig (se_mat_in,xlabels,ylabels):
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  ax.set_aspect('equal')
  plt.imshow(np.array(se_mat_in),cmap=plt.cm.Blues,interpolation='nearest',origin='lower') #interpolation='nearest')#, cmap=plt.cm.ocean)
  plt.grid(True)
  plt.xlabel("Side effects")
  plt.ylabel("Durgs")
  ax.set_xticks=[range(len(xlabels))]
  ax.set_yticks=[range(len(ylabels))]
  ax.set_xticklabels(xlabels, rotation=45, rotation_mode="anchor")
  ax.set_yticklabels(ylabels, rotation=45, rotation_mode="anchor")
  plt.savefig('app/static/XXX_DrugSideEffects.png', dpi=150)
  #plt.show()
#############################
def PlotBarh (top_se_names, top_se_prob,rand_i):
  plt.rcParams.update({'font.size': 18})
  y_pos = np.arange(len(top_se_names))
  plt.clf()
  plt.cla()
  plt.autoscale()
  plt.barh(y_pos, top_se_prob)#,  alpha=0.4) #
  plt.yticks(y_pos, top_se_names)
  plt.xlabel('Probability')
  plt.title('Top Predicted Side Effects')

  #plt.savefig('static/XXX_DrugSideEffects'+rand_i+'.png', dpi=150,bbox_inches='tight')
  plt.savefig('static/XXX_DrugSideEffects.png', dpi=150,bbox_inches='tight')
  plt.close
  #plt.show()


#############################
def PlotVenn (s):
  plt.clf()
  plt.cla()
  plt.autoscale()
  plt.rcParams.update({'font.size': 22})
  plt.title('Predicted and Known Side Effects')

  # Set values 
  v = venn2(subsets=s, set_labels=('Known', 'Predicted'))
  # Subset labels
  v.get_label_by_id('10').set_text(str(s['10']))
  v.get_label_by_id('01').set_text(str(s['01']))
  #v.get_label_by_id('11').set_text(str(s['11']))

  # Subset colors
  v.get_patch_by_id('10').set_color('green')
  v.get_patch_by_id('01').set_color('salmon')
  #v.get_patch_by_id('11').set_color('blue')

  # Subset alphas
  v.get_patch_by_id('10').set_alpha(0.4)
  v.get_patch_by_id('01').set_alpha(1.0)
  
  if (s['11'] > 0 ): # if no overlap 
      v.get_label_by_id('11').set_text(str(s['11']))
      v.get_patch_by_id('11').set_color('blue')
      v.get_patch_by_id('11').set_alpha(0.7)

  # Border styles
  #c = venn2_circles(subsets=s, linestyle='solid')
  #c[0].set_ls('dashed')  # Line style
  #c[0].set_lw(2.0)       # Line width
  #plt.show()
  plt.savefig('static/XXX_VennDiagram.png', dpi=150,bbox_inches='tight')
#############################
@app.before_request
def before_request():
    
    if not ('BpredSE' in globals() ) :
      global BpredSE, Bscore
      BpredSE = pk.load(open(pred_file,'rb'))
      Bscore  = pk.load(open(score_file,'rb'))

####

@app.after_request
def add_header(response):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    response.headers['X-UA-Compatible'] = 'IE=Edge,chrome=1'
    response.headers['Cache-Control'] = 'public, max-age=0'
    return response


@app.route('/')
@app.route('/index')
@app.route('/input')
#def index():
#		user={'nickname':'Aqeel'}
#		return render_template("index.html",title="Home",user=user)
def drug_input():
  # Get drug_list
  drug_list = GetDrugList ()
  return render_template("input.html",drug_list=drug_list)

@app.route("/data")
@app.route("/data/<drug_name>")
## @app.route("/data/<search_input>")
def data(drug_name):
    set_TOP = 15
    # get drug_id and its data 
    #drug_id   = request.args.get('drug_id')
    global BpredSE, Bscore
    #print "=====++++",drug_id
    #drug_name = GetDrugName (drug_id)[0]['drug_name']
    #print "====++++",drug_name
    TOP_predicted_se = BpredSE['predSE'][drug_name]['predicted_se'][:set_TOP]
    TOP_Y_probability= BpredSE['predSE'][drug_name]['Y_probability'][:set_TOP]
    #print TOP_predicted_se
    #global TOP_predicted_se,TOP_Y_probability
    #http://nvd3.org/examples/multiBarHorizontal.html
    JsonProb =   {
        "key": "Predicted probability of side effect",
        "color": "#d62728"#,
        #"values": [{'label':'se1', 'value':22},{'label':'se2','value':44}]
        }
    
    values = []
    #print len(TOP_predicted_se)
    for i in range (len(TOP_predicted_se)):
        val = {}
        val["label"] = TOP_predicted_se[i]
        val["value"] = TOP_Y_probability[i]
        values.append(val)
    JsonProb['values'] = values
    #print json.dumps([JsonProb])
    return json.dumps([JsonProb])

@app.route('/output')
def LogRegPredict ():
  drug_id   = request.args.get('drug_id')
  drug_name = GetDrugName (drug_id)[0]['drug_name']
  #rand_i = str(time.time())
  #print drug_name
  print "++++",drug_id

  # Get drug_list
  drug_list = GetDrugList ()
  #BpredSE = joblib.load(pred_file)
  #Ascore  = joblib.load(score_file)

  #BpredSE = pk.load(open(pred_file,'rb'))
  #Bscore  = pk.load(open(score_file,'rb'))
  global BpredSE, Bscore
  predSE      = BpredSE['predSE']
  drug_names  = BpredSE['drug_names']
  se_names    = BpredSE['se_names']
  Ascore      = Bscore['Ascore']
  #
  #print Ascore

  #print 'BpredSE:',BpredSE se_names_top Y_pred, Y_probability
  # Give error if drug not not found 
  #print 'drug_names',drug_names
  if (drug_name not in predSE):
    return ('</br> </br> <b>Sorry, the drug is not present in our database!</b>')

  # Get data for the drug 
  predicted_se  = predSE[drug_name]['predicted_se']  # this is sorted 
  Y_pred        = predSE[drug_name]['Y_pred']        # this contains all (not sorted)
  Y_probability = predSE[drug_name]['Y_probability'] # this is sorted
  predicted_se_count  = len(predicted_se)
  print 'predicted_se_count',predicted_se_count
  # Get Roc_AUC
  roc_auc = []
  for see in predicted_se : 
    roc_auc.append (Ascore[see]['roc_auc'])


  ###########
  # sort by Y_probability and get >.5 
  #Y_prob_ind =  sorted(range(len(Y_probability)),key=lambda x:Y_probability[x],reverse=True)[:predicted_se_count]
  #top_se_names = [se_names[i] for i in Y_prob_ind[:20]]
  #other_se_names = [se_names[i] for i in Y_prob_ind[:50]]
  #top_se_prob  = [Y_probability[i] for i in Y_prob_ind[:20]  ]  #[:20]
  #other_se_prob  = [Y_probability[i] for i in Y_prob_ind[:50]  ]
  #############
  #print [Y_probability[i] for i in Y_prob_ind  ]
  #print Y_prob_ind
  #print top_se_names

  #global TOP_predicted_se   
  #TOP_predicted_se = predicted_se[:15]
  #global TOP_Y_probability  
  #TOP_Y_probability = Y_probability[:15]

  #PlotBarh (predicted_se[:20], Y_probability[:20],rand_i)
  #PlotBarh (predicted_se[:20], roc_auc[:20],rand_i)


  #accuracy      = Ascore [drug_name]['accuracy']
  #roc_auc       = format(Ascore [drug_name]['roc_auc'], '.2f') 
  #top_se_count  = len(top_se)

  # get predicted se for the drug 
  #predicted_se        = [top_se[i] for i in range(top_se_count) if Y_pred[i] == 1]

  print "len:predicted_se",predicted_se_count#, #predicted_se
  if (predicted_se_count ==  0):
    return ('</br> </br> <b>Due to limited dataset, no side effect was predicted! Please go back and try another drug.</b>')

  # get known se for the drug
  known_se = GetSEforForTheDrug(drug_id)
  known_se_count = len(known_se)
  pred_se_known = ['Yes' if see in known_se else 'No' for see in predicted_se]
  pred_se_known_count =  sum([1 if see in known_se else 0 for see in predicted_se])
  # Calculate overlap/match between top and known 
  #se_initial_count = sum([1 if se in top_se else 0 for se in known_se])
  #print "se_initial_count:",se_initial_count

  # Calculate overlap/match between predicted and known 
  #se_match_count    = sum([1 if se in predicted_se else 0 for se in known_se])
  #se_accuracy_count = sum([1 if i in range(predicted_se_count) else 0 for i in range(top_se_count)])
  #print "se_match_count:",se_match_count
  # 
  #c = {'10':int(known_se_count)-int(pred_se_known_count),'01':int(predicted_se_count)-int(pred_se_known_count),'11':0}
  c = {'10':int(known_se_count)-int(pred_se_known_count),'01':int(predicted_se_count)-int(pred_se_known_count),'11':int(pred_se_known_count)}
  PlotVenn  (c)
  # 

  return render_template("output.html",drug_list =drug_list, drug_name=drug_name, predicted_se=predicted_se,
   pred_se_prob=Y_probability, roc_auc = roc_auc,#) #,
    pred_se_known=pred_se_known, pred_se_known_count=pred_se_known_count, known_se_count=known_se_count)
      #se_match_count=se_match_count,se_initial_count=se_initial_count, known_se_count=known_se_count, top_se_count=top_se_count, 
      #predicted_se_count=predicted_se_count,accuracy=accuracy,roc_auc=roc_auc)


if __name__ == '__main__' :
  app.run(host='0.0.0.0',port=5000,debug=True)