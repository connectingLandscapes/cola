# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 21:00:16 2023
Connects EE from command line
@author: ig299
"""

#%%
# IMPORTS
import ee
import csv
import sys 
import pandas as pd


def list_assets(parent_asset):
    """
    List all assets in a specified Earth Engine folder or project.

    Parameters
    ----------
    parent_asset : str
        Example:
        - "projects/my-project/assets"
        - "users/username"
        - "projects/earthengine-legacy/assets/users/username"
    """

    try:
        result = ee.data.listAssets({"parent": parent_asset})
        assets = result.get("assets", [])

        if not assets:
            print("No assets found.")
            return

        #print(f"\nAssets under: {parent_asset}\n")

        #for asset in assets:
        #    print(f"Name : {asset['name']}")
        #    print(f"Type : {asset['type']}")
        #    print("-" * 60)

    except Exception as e:
        print(f"Error: {e}")


# ---------------------------------------------------------------------
# Recursively list assets
# ---------------------------------------------------------------------
def walk_assets(parent, assets_list, level=0):
    """
    Recursively walks through an Earth Engine asset folder.

    Parameters
    ----------
    parent : str
        Parent asset path.
    assets_list : list
        List where asset information will be stored.
    level : int
        Recursion depth.
    """

    try:
        result = ee.data.listAssets({"parent": parent})
        for asset in result.get("assets", []):
            asset_info = {
                "Level": level,
                "Type": asset["type"],
                "Name": asset["name"],
                "Parent": parent
            }
            assets_list.append(asset_info)
            #print(
            #    f"{'    '*level}"
            #    f"{asset['type']:<18} {asset['name']}"
            #)
            # Continue walking folders and image collections
            if asset["type"] in ["FOLDER", "IMAGE_COLLECTION"]:
                walk_assets(asset["name"], assets_list, level + 1)
                #
            #
    except Exception as e:
        print(f"Could not access {parent}: {e}")

# ---------------------------------------------------------------------
# Export to CSV
# ---------------------------------------------------------------------
def export_csv(asset_list, output_file):
    with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(
            csvfile,
            fieldnames=["Level", "Type", "Parent", "Name"]
        )
        writer.writeheader()
        for row in asset_list:
            writer.writerow(row)
    print(f"\nInventory exported to:\n{output_file}")

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------


def main() -> None:
    #%%
    # Start timer
    #
    # INPUTS
    # Path to file holding xy coordinates
    param1 = sys.argv[1] # project name # param1 = 'gonzalezivan'
    param2 = sys.argv[2] # tempfolder were to write # param2 = 'C:/cola/colaNAR2026030412391405/'
    #
    # C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_connectEE.py" gonzalezivan C:/cola/colaFLF2026072014371305/
    # param1 = 'gonzalezivan'
    # param2 = 'C:/cola/cola2/id2000.txt'
    # Convert distance threshold to float or integer
    try:
        ee.Authenticate()
        ee.Initialize(project = param1)
        print('\n Cola2: EE initialized')
    except:
        print('\n Cola2 ERROR: no EE initialized\n')
        exit(1)   
    #
    # Convert gaussian smoother radius to integer
    try:
        collection = ee.ImageCollection('NOAA/DMSP-OLS/NIGHTTIME_LIGHTS')
        collection_id = collection.getMapId().get('mapid')
        print('\n Cola2: EE conection works!\n')
        with open(param2+'/ee_conected.txt', "w") as file:
          file.write(collection_id )
          # print(collection.getMapId())
        #
        tasks = ee.data.listOperations()
        rows = [''] * len(tasks)
        counter = 0
        for task in tasks:
            #print(f"Name: {task.get('name')}")
            #print(f"Metadata: {task.get('metadata')}")
            #print("-" * 40)
            md = task.get('metadata')
            ec = md.get("batchEecuUsageSeconds") or md.get("eecuUsageSeconds") or md.get("batchEecuSeconds") or "0"
            dd = {'name': task.get('name', ''),
                  'typeMetadata': md.get('@type'),
                    'state': md.get('state', ''),
                    'description': md.get('description', ''),
                    'priority': md.get('priority', ''),
                    'createTime': md.get('createTime', ''),
                    'updateTime': md.get('updateTime', ''),
                    'startTime': md.get('startTime', ''),
                    'endTime': md.get('endTime', ''),
                    'typeTask': md.get('type', ''),
                    'destinationUris': md.get('destinationUris', ' ')[0],
                    'attempt': md.get('attempt', ''),
                    'progress': md.get('progress', ''),
        # =============================================================================
        #               'displayName': md.get('stages', '')[0].get('displayName', ''),
        #               'completeWorkUnits': md.get('stages')[0].get('completeWorkUnits', ''),
        #               'totalWorkUnits': md.get('stages')[0].get('totalWorkUnits', ''),
        #               'description': md.get('stages')[0].get('description', ''),
        #               'displayName2': md.get('stages')[1].get('displayName', ''),
        #               'completeWorkUnits2': md.get('stages')[1].get('completeWorkUnits', ''),
        #               'totalWorkUnits2': md.get('stages')[1].get('totalWorkUnits2', ''),
        #               'description2': md.get('stages')[1].get('description2', ''),
        # =============================================================================
                'eecu': float( ec ),
                'done': task.get('done', ''),
                'response': task.get('response', ''),
                'error': task.get('error', '')
                }
            rows[counter] = dd
            counter += 1
        # -----------------------------------------------------------------------------
        # Create DataFrame
        # -----------------------------------------------------------------------------
        df = pd.DataFrame(rows)
        # Convert timestamps to datetime
        df["createTime"] = pd.to_datetime(df["createTime"], errors="coerce")
        df["updateTime"] = pd.to_datetime(df["updateTime"], errors="coerce")
        df["startTime"] = pd.to_datetime(df["startTime"], errors="coerce")
        df["endTime"] = pd.to_datetime(df["endTime"], errors="coerce")
        # -----------------------------------------------------------------------------
        # Export all tasks to CSV
        # -----------------------------------------------------------------------------
        pathh = param2+ "/gee_tasks_eecu.csv"
        df.to_csv(pathh, index=False)
        print(f"Tasks saved: {pathh}")
        # -----------------------------------------------------------------------------
        # Monthly EECU summary
        # -----------------------------------------------------------------------------
        monthly = (
            df.assign(Month=df["endTime"].dt.to_period("M"))
              .groupby("Month", as_index=False)["eecu"]
              .sum()        )
        monthly["Month"] = monthly["Month"].astype(str)
        monthly.to_csv(param2 + "/gee_eecu_monthly_summary.csv", index=False)
        print("\nMonthly EECU Consumption")
        print(monthly)
        # List assets
        ROOT_ASSET = "projects/" + param1 + "/assets"
        assets = []
        walk_assets(ROOT_ASSET, assets)
        path2 = param2 + "/gee_assets_inventory.csv"
        export_csv(assets, path2)
        print(f"\nTotal assets found: {len(assets)}")
        print(f"Assets saved: {pathh}")
        #
    except ValueError as e:
        print('Cola 2 ERROR: ', e)
        exit(1)


if __name__ == "__main__":
    main()
