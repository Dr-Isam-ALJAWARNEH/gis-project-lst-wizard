# 🌍 LST Tool for ArcGIS Pro

A Python-based ArcGIS Pro script tool for calculating **Land Surface Temperature (LST)** from **Landsat 7, Landsat 8, and Landsat 9** satellite imagery. The tool automates NDVI, emissivity, brightness temperature, and final LST extraction, with support for both **single scene processing** and **batch mode**.

## 🚀 Features

- Supports Landsat 7, 8, and 9
- Detects metadata and sensor type automatically using `MTL.txt`
- Outputs:
  - Land Surface Temperature (°C)
  - NDVI *(optional)*
  - Emissivity *(optional)*
  - Brightness Temperature *(optional)*
- Supports single scene or batch mode
- Optional AI summary using OpenAI API
- Skips failed scenes without stopping

📦 Repository Contents

| File/Folder | Description |
|------------|-------------|
| `LST_Tool.py` | Main Python script |
| `LST.atbx` | ArcGIS Pro toolbox |
| `README.md` | Documentation |
| Sample Output| Example output files |

🔧 Requirements

Software
- ArcGIS Pro (Python 3 environment)
- Spatial Analyst Extension

Python Libraries
| Module | Required | Included |
|--------|----------|----------|
| `arcpy` | ✔ | Yes |
| `arcpy.sa` | ✔ | Yes |
| `requests` | Optional | No |


📂 Input Requirements (Landsat 7, 8, or 9 Collection 2 Level 1 from USGS) study area of your choice.
https://earthexplorer.usgs.gov/
Or you can use data available in the following the link.
https://drive.google.com/drive/folders/1gAoYpVjcApCFRLJRA9iU8mRX-RfzXPak?usp=sharing

Expected Landsat scene format:
SCENE_FOLDER/
│── *_MTL.txt 
│── Bands (TIF format)

 ▶ Running in ArcGIS Pro

1. Open ArcGIS Pro  
2. Go to **Catalog → Add Toolbox**  
3. Select `LST.atbx`  
4. Run the tool and choose settings  

