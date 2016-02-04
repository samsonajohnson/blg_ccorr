import dill
import numpy as np
import scipy
import scipy.interpolate
import ipdb
import pickle

class test_class:
    def __init__(self,arg1,arg2,arg3,filename = './test.pickle'):
        self.arg1 = arg1
        self.arg2 = arg2
        self.arg3 = arg3
        self.x = np.linspace(0,2*np.pi,100)
        self.y = np.sin(self.x)
        self.testspline = scipy.interpolate.interp1d(self.x,self.y,kind='cubic')

    def foo(self):
        print self.arg1
    def bar(self):
        print self.arg1+self.arg2

    def pickle(self,filename = './test.pickle'):
        with file(filename,'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

    def unpickle(self,filename = './test.pkl'):
       with file('test,pkl', 'rb') as f:
            return pickle.load(f)
        


                     
