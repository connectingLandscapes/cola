# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:08:33 2026

@author: gonza
"""

import ee
import csv

# ---------------------------------------------------------------------
# Authenticate and initialize Earth Engine
# ---------------------------------------------------------------------
def initialize_gee():
    try:
        ee.Initialize()
        print("Earth Engine initialized.")
    except Exception:
        ee.Authenticate()
        ee.Initialize()
        print("Earth Engine authenticated and initialized.")



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

        print(f"\nAssets under: {parent_asset}\n")

        for asset in assets:
            print(f"Name : {asset['name']}")
            print(f"Type : {asset['type']}")
            print("-" * 60)

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
            print(
                f"{'    '*level}"
                f"{asset['type']:<18} {asset['name']}"
            )
            # Continue walking folders and image collections
            if asset["type"] in ["FOLDER", "IMAGE_COLLECTION"]:
                walk_assets(asset["name"], assets_list, level + 1)

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

if __name__ == "__main__":

    initialize_gee()
    ee.Initialize(project = 'gonzalezivan')
    # Change this to your asset folder
    ROOT_ASSET = "projects/gonzalezivan/assets"
    parent_asset = "projects/gonzalezivan/assets"
    parent = "projects/gonzalezivan/assets"

    # list_assets(parent_asset)
    # Replace with your asset root
    assets = []
    walk_assets(ROOT_ASSET, assets)
    assets[0]
    export_csv(assets, "C:/cola/gee_assets_inventory.csv")
    print(f"\nTotal assets found: {len(assets)}")

