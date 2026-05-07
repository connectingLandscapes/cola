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


def checkMyTasksSeconds ():
    mytasks = ee.batch.Task.list()
    eecu_secs = 0
    for e1 in mytasks:
        eecu_secs  += e1.status().get('batch_eecu_usage_seconds')
        # print( eecu_secs  )
    #    
    mytasks2 = ee.data.listOperations()
    eecu_secs2 = 0
    for e2 in mytasks2:
        eecu_secs2 += float(e2.get('metadata').get('batchEecuUsageSeconds'))
        #print( eecu_secs2 )
    #
    return eecu_secs2
#
#
def asset_exists(asset_id):
    try:
        ee.data.getAsset(asset_id)
        return True
    except ee.EEException:
        return False
