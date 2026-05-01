# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 09:46:37 2026

@author: ig299
"""

def taskTodict(mytask):
    conf = mytask.config
    #    
    taskdict = {
        "id0" : str(mytask.id),
        "name": str(mytask.name),
        "workload_tag": str(mytask.workload_tag),
        "task_type": str(mytask.task_type),
        "state": str(mytask.state),
        #"expressionv": str(conf.get('expression')),
        "description": str(conf.get('description')),
        "destination" : str(conf.get('assetExportOptions').get('earthEngineDestination').get('name')),
        "overwrit": str(conf.get('assetExportOptions').get('earthEngineDestination').get('overwrite'))
        }
    return taskdict

def writeTask(task, folder):
    # folder = 'C:/Users/ig299/cola/eelogpath'
    import pandas as pd
    taskdict = taskTodict(task)
    outFile = folder + "/eelog_" + taskdict.get('id0') + '.csv'    
    df = pd.DataFrame(taskdict , index = [0])
    # pd.DataFrame([taskdict])
    df.to_csv(outFile, index=False)


def checkMyTasks ():
    mytasks = ee.batch.Task.list()
    for t in mytasks:
        print(t.status())
    mytasks2 = ee.data.listOperations()
    return 1
    #len(mytasks2)
#