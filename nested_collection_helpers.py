import random
import IPython
from collections import defaultdict
import pprint as pp

# another appoach is to NOT use nestedDitc! 
# instead use a single dict (maybe a default)
# but with tuple keys.  The keys can the have a hierarchy
# mykey = 'Bearing', 'P', 4, 'G4' 



class NestedDict(dict):    
    """Use of __missing__ will create keys from any typo"""
    def __missing__(self, key):  
        self[key] = NestedDict()
        return self[key]

    
def walkDict_values( inDict, depth=0, retDict={}):
    retStr = ""
    """walk all the NestedDict and dict"""
    pre = "\t"*depth
    for k, v in inDict.items():
        if not isinstance(v, (NestedDict, dict)):
            retStr += (f"{pre}{k} -> {type(v)}\n")
            if type(v) not in retDict:
                retDict[type(v)] = list()
            retDict[type(v)].append(v)
        else:
            print(f"{pre}{k}\n")
            walkDict_values(v, depth=depth +1, retDict=retDict )
    depth -= 1
    return retDict


def find_key(d:dict, value):
    """walks a nested dict looking for a value. Returns an order list of
    of keys that references the value"""
    for k,v in d.items():
        if isinstance(v, dict):
            p = find_key(v, value)
            if p:
                return [k] + p
        elif v == value:
            return [k]


def rget(obj, attrstr, default=None, delim='.'):
    """Recursive get of nested containers, dict, list, tuple object covered.
    call this:
    value = rget(someobj, 'foo.bar.baz')"""
    try:
        parts = attrstr.split(delim, 1)
        attr = parts[0]
        attrstr = parts[1] if len(parts) == 2 else None
        if isinstance(obj, dict): value = obj[attr]
        elif isinstance(obj, list): value = obj[int(attr)]
        elif isinstance(obj, tuple): value = obj[int(attr)]
        elif isinstance(obj, object): value = getattr(obj, attr)
        if attrstr: return rget(value, attrstr, default, delim)
        return value
    except Exception:
        return default




class DefaultBunch(defaultdict):
    """Bunch and munch are classes that have been developed to allow .atrribute indexing.
    A Bunch is a dictionary-like object that allows you to access its keys as attributes. 
    This means you can use dot notation (e.g., my_bunch.key) instead of the traditional 
    dictionary access (e.g., my_bunch['key']).
    A Bunch is a dictionary that provides attribute-style access (a la JavaScript objects).
    see Scikit-learn
    
    Useage:
    flexbunch = lambda: DefaultBunch(flexbunch)
    f = flexbunch()
    f.easily.create.deeply.nested.structures.using.dot.notation = 1
    """
   
    def __getattr__(self, item): return self.__getitem__(item)
    """If the attribute exists __getitem__ it"""

    def __setattr__(self, item, value): return self.__setitem__(item, value)
    """return self.__setitem__(item, value)"""



if __name__ == "__main__":




    t = NestedDict()
    c='curves'
    s='segments'
    Ch='chains'
    Ra = 'rays'
    p='pts'
    t[s]= []  #list of segments
    pts = []
#set CL
    t["CL-K96"] = random.random()
#make Chain of road
    t[Ch]["Road"] = random.random()
    ref = t[Ch]["Road"]
#set ProGr
    t[Ch]["ProGr"] =  random.random()
#t[Ch]["GL_A"]
    for myPG in range(1,6):
        t["ProGr"][myPG] = random.random()


    produceData = True 

    if produceData:
        t["GL"]["A"] = random.random()
        t["GL"]["P_1"]["G_A"] = random.random()
        t["GL"]["P_5"]["G_A"] = random.random()
        t["GL"]["B"] = random.random()
        t["GL"]["C"] = random.random()
        t["GL"]["D"] = random.random()
        t["Bearing"]["P_1"]["G_D"] = random.random()
        t["Bearing"]["P_2"]["G_D"] = random.random()
        t["Bearing"]["P_3"]["G_D"] = random.random()
    for GL in ["A", "B", "C"]:
        for i in range(1,6):
            #set many of the GL Pier Points
            t["Bearing"]["P_" + str(i)]["G_" + GL] = random.random()

    if produceData:
        t["Bearing"]["P_5"]["G_C"] = random.random()
        t["Bearing"]["P_4"]["G_C"] = random.random()
        t["Bearing"]["P_3"]["G_C"] = random.random()
        t["Bearing"]["P_1"]["G_G"] = random.random()
        t["Bearing"]["P_2"]["G_G"] = random.random()
        t["Bearing"]["P_3"]["G_G"] = random.random()
        t["Bearing"]["P_4"]["G_G"] = random.random()
        t["Bearing"]["P_5"]["G_G"] = random.random()
        t["Bearing"]["P_1"]["G_J"] = random.random()
        t["Bearing"]["P_2"]["G_J"] = random.random()
        t["Bearing"]["P_3"]["G_H"] = random.random()
        t["Bearing"]["P_3"]["G_J"] = random.random()
        t["Bearing"]["P_4"]["G_H"] = random.random()
        t["Bearing"]["P_4"]["G_J"] = random.random()
        t["Bearing"]["P_5"]["G_H"] = random.random()
        t["Bearing"]["P_5"]["G_J"] = random.random()
        
#print( t.walkDict() )

    print( f"{len(walkDict_values(t))=}" )
    for k,mylist in walkDict_values(t).items():
        for obj in mylist:
            print( obj.__repr__(), find_key(t, obj))

    flexbunch = lambda: DefaultBunch(flexbunch)
    f = flexbunch()
    f.easily.create.deeply.nested.structures.using.dot.notation = 1
    print("\n",f)

    import IPython
    print("\n\nNow start IPython.embed")
    IPython.embed()
