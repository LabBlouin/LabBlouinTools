# ------------------------------ #
# svm.py                         #
# ------------------------------ #
# Winter 2014                    #
# ------------------------------ #
# Alex Safatli                   #
# ------------------------------ #
# SVM Managing, Cross-Validation #
# ------------------------------ #

import svmutil, random, os
from numpy import arange
#from collections import Counter as C 

''' SVM Model                  '''

class svm_model():

    def __init__(self, data=[], labels=[], kernel=svmutil.RBF, c=10):
        self.__svmparam__             = svmutil.svm_parameter('-q')
        self.__svmparam__.kernel_type = kernel
        self.__svmparam__.C           = c
        self.c                        = c
        self.data                     = data
        self.labels                   = labels
        if len(data) > 0:
            self.problem = svmutil.svm_problem(labels,data)
            self.model   = svmutil.svm_train(self.problem,self.__svmparam__,'-q')
        else:
            self.problem = None
            self.model   = None

    def optimize(self,low=1,up=50,steps=5,by='proportions'):
        '''
        Optimize c-value on the basis of by (default:proportions of correct classifications).
        '''
        optim = {}
        print "Optimizing the C-value:"
        for c_val in arange(low,up,steps):
            model = svm_model(self.data,self.labels,c=c_val)
            valid = cross_validate(model)
            valid.perform()
            if by == 'proportions': optim[c_val] = valid.proportion
            else:                   optim[c_val] = valid.meanFscore

        max_c = max(optim,key=lambda d: optim[d])

        self.c              = max_c
        self.__svmparam__.C = max_c
        self.model = svmutil.svm_train(self.problem,self.__svmparam__,'-q')
        print "Optimum C: %f"%(max_c)

    def predict(self,data):

        # Predict.
        predict,_,_ = svmutil.svm_predict([-1]*len(data),data,self.model)
        return predict

    def load(self,fin):
        self.model = svmutil.svm_load_model(fin)

    def save(self,fin):
        o = open(fin,'w')
        o.close()        
        svmutil.svm_save_model(fin,self.model)


''' Cross-Validation           '''

class cross_validate():

    def __init__(self, model, test=[]):
        self.model  = model
        self.data   = model.data
        self.labels = model.labels
        self.test   = test

    def _randomize(self):
        '''
        Randomize all data.
        '''
        indicies = range(len(self.labels))
        random.shuffle(indicies)
        self.data_r =   [self.data[x]   for x in indicies]
        self.labels_r = [self.labels[x] for x in indicies]

    def _split(self):
        '''
        Splits all data into n
        lists.
        '''
        n = int(len(self.data)*0.1)
        if n == 0: n = 1        
        self.splits = []
        l1 = self.data_r
        l2 = self.labels_r
        for i in xrange(0,len(self.data),n):
            self.splits.append((l1[i:i+n],l2[i:i+n]))
        self.full = []
        for i in xrange(len(self.data)):
            self.full.append((l1[i],l2[i]))        

    def _exclude(self, split_list):
        '''
        Returns a list of all items
        but those found in a split
        list.
        '''
        check_list = []
        for item in split_list[0]: check_list.append(item)
        exclude_out = []
        for item in self.full:
            if item[0] not in check_list: exclude_out.append(item)
            else: check_list.remove(item[0])
        return exclude_out

    def _test(self, split, exclude):
        '''
        Tests all items in a split
        part of the data across
        everything else. Returns 
        all relevant statistics
        for F-score evaluation.
        '''
        falseNeg = 0
        falsePos = 0
        truePos = 0
        excludeLabels = []
        excludeData   = []
        for item in exclude:
            excludeLabels.append(item[1])
            excludeData.append(item[0])
        svm_m     = svm_model(excludeData,excludeLabels,
                              kernel=self.model.__svmparam__.kernel_type,
                              c=self.model.c)
        split_len = len(split[0])
        #for i in xrange(split_len):
        predict,acc,val = svmutil.svm_predict(split[1],split[0],svm_m.model)
        #print 'Predicted label: %d\t Real label:%d'%(int(predict[0]),int(split[1][i]))
        prediction, reality = predict, split[1]
        return (prediction, reality)
    
    def _confusionMatrixCounter(self, stats):
        ''' compute the counts for each label'''
        matrix={}
        classes = list(set(self.labels))
        for e in classes:
            matrix[e]={}
            for f in classes:
                matrix[e][f] = 0
        for s in stats:
            reality, prediction = s[1], s[0]
            for l in set(reality):
                for i in range(len(reality)):
                    if reality[i] != l:
                        continue
                    else:
                        matrix[l][prediction[i]] += 1
        return matrix

    def _processStats(self, stats):
        '''
        Process a stats list. Get F-score.
        '''
        cvMatrix = self._confusionMatrixCounter( stats )
        '''cvMatrix = {}
        for e in stats[0].keys():
            cvMatrix[e]={}
            for f in stats[0].keys():
                cvMatrix[e][f]=0
        for item in stats:
            for k,v in item.iteritems():
                for cls, val in v.iteritems():
                    cvMatrix[k][cls]+=val'''
        # specificity Sp = TP / (TP+FP)
        # sensitivity Sv = TP / (TP+FN)
        # F score = 2 * (Sp*Sv)/(Sp+Sv)        
        self.f1stats={}
        meanFscore=0
        correctlyclass=0
        total=0
        for i in cvMatrix.iterkeys():
            falseNeg = 0
            falsePos = 0
            truePos  = 0
            trueNeg  = 0            
            self.f1stats[i] = []
            for k, v in cvMatrix.iteritems():
                for ke, va in v.iteritems():
                    if i == k == ke: 
                        truePos  = cvMatrix[k][ke]
                        correctlyclass += cvMatrix[k][ke]
                    elif i == k and k != ke: falseNeg += cvMatrix[k][ke]
                    elif i == ke and k != i: falsePos += cvMatrix[k][ke]
                    else: trueNeg += cvMatrix[k][ke]
            #specificity
            sp=(float(truePos)+1) / ((truePos+falsePos) +1)
            self.f1stats[i].append(sp)
            #sensitivity
            sn= (float(truePos)+1) / ((truePos+falseNeg)+1)
            self.f1stats[i].append(sn)  
            #F-score
            fscore = (2*(sp*sn)+1)/((sp+sn)+1)
            self.f1stats[i].append(fscore)
            meanFscore += fscore
            #correctlyclass+= truePos+trueNeg
            total+= truePos+falsePos+trueNeg+falsePos
        self.meanFscore = meanFscore/len(self.f1stats)            
        self.stats = cvMatrix
        self.proportion = float(correctlyclass)/len(self.data)
        


    def perform(self):
        '''
        Perform the cross-validation.
        '''
        line = "10% Cross-validation with optimize model"
        print "#"*(len(line)+2)+'\n'+'#'+line+'#'+'\n'+"#"*(len(line)+2)        
        stats = []
        self._randomize()
        self._split()
        for s in self.splits:
            # s is a chunk of the full set of data.
            e = self._exclude(s) # List of everything else.
            stats.append(self._test(s,e))
        # Combine statistics and get F-score.
        self._processStats(stats) 