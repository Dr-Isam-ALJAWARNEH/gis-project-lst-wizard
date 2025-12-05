import arcpy
import os
import datetime
import textwrap
from arcpy.sa import *  # Spatial Analyst

try:
    import requests
except ImportError:
    requests = None

# Messaging Helpers
def log_message(msg):
    try: arcpy.AddMessage(msg)
    except: print(msg)

def log_warning(msg):
    try: arcpy.AddWarning(msg)
    except: print("WARNING:", msg)

def log_error(msg):
    try: arcpy.AddError(msg)
    except: print("ERROR:", msg)


# ----------------------------------------------------------
# API Key Handler
# ----------------------------------------------------------
def get_api_key(tool_param_key):
    if tool_param_key:
        return tool_param_key

    return os.environ.get("OPENAI_API_KEY") or os.environ.get("LST_TOOL_OPENAI_KEY")


# ----------------------------------------------------------
# Landsat Metadata Helpers
# ----------------------------------------------------------
def find_mtl_file(folder):
    for f in os.listdir(folder):
        if f.upper().endswith("MTL.TXT"):
            return os.path.join(folder, f)
    return None


def parse_mtl(path):
    m = {}
    with open(path, "r") as f:
        for line in f:
            if "=" in line:
                key, val = line.split("=", 1)
                m[key.strip()] = val.strip().strip('"')
    return m


def detect_landsat_sensor(meta):
    s = meta.get("SPACECRAFT_ID", "").upper()
    if "LANDSAT_7" in s: return "L7"
    if "LANDSAT_8" in s: return "L8"
    if "LANDSAT_9" in s: return "L9"

    scene = meta.get("LANDSAT_SCENE_ID", "").upper()
    if scene.startswith(("LE07", "LT07")): return "L7"
    if scene.startswith(("LC08", "LO08")): return "L8"
    if scene.startswith(("LC09",)): return "L9"
    return None


def find_band_file(folder, suffix):
    suffix = suffix.upper()
    for f in os.listdir(folder):
        if f.upper().endswith(suffix):
            return os.path.join(folder, f)
    return None


# ----------------------------------------------------------
# Main LST Calculation
# ----------------------------------------------------------
def compute_landsat_lst_for_scene(scene_folder, thermal_band_number=None,
                                  out_lst_path=None, save_ndvi=False,
                                  save_emissivity=False, save_bt=False):

    log_message(f"Processing scene folder: {scene_folder}")

    mtl_path = find_mtl_file(scene_folder)
    if not mtl_path:
        raise arcpy.ExecuteError(f"No MTL file found in {scene_folder}")

    mtl = parse_mtl(mtl_path)
    sensor = detect_landsat_sensor(mtl)
    if not sensor: raise arcpy.ExecuteError("Unsupported Landsat mission.")

    log_message(f"Detected sensor: {sensor}")

    # ---- Band Assignment ----
    if sensor in ("L8", "L9"):

        tb = thermal_band_number if thermal_band_number in (10, 11) else 10
        log_message(f"Using thermal band {tb} for Landsat {sensor}.")

        red_path = find_band_file(scene_folder, "_B4.TIF")
        nir_path = find_band_file(scene_folder, "_B5.TIF")
        thr_path = find_band_file(scene_folder, f"_B{tb}.TIF")

        rad_mult = float(mtl[f"RADIANCE_MULT_BAND_{tb}"])
        rad_add = float(mtl[f"RADIANCE_ADD_BAND_{tb}"])
        k1 = float(mtl[f"K1_CONSTANT_BAND_{tb}"])
        k2 = float(mtl[f"K2_CONSTANT_BAND_{tb}"])
        lambda_thermal = 10.895e-6 if tb == 10 else 12.0e-6

        refl_mult_red = float(mtl.get("REFLECTANCE_MULT_BAND_4", "1"))
        refl_add_red = float(mtl.get("REFLECTANCE_ADD_BAND_4", "0"))
        refl_mult_nir = float(mtl.get("REFLECTANCE_MULT_BAND_5", "1"))
        refl_add_nir = float(mtl.get("REFLECTANCE_ADD_BAND_5", "0"))

    else:  # Landsat 7
        log_message("Using thermal band 6 for Landsat 7.")

        red_path = find_band_file(scene_folder, "_B3.TIF")
        nir_path = find_band_file(scene_folder, "_B4.TIF")
        thr_path = find_band_file(scene_folder, "_B6_VCID_1.TIF") or find_band_file(scene_folder, "_B6.TIF")

        rad_mult = float(mtl["RADIANCE_MULT_BAND_6_VCID_1"])
        rad_add = float(mtl["RADIANCE_ADD_BAND_6_VCID_1"])
        k1 = float(mtl["K1_CONSTANT_BAND_6_VCID_1"])
        k2 = float(mtl["K2_CONSTANT_BAND_6_VCID_1"])
        lambda_thermal = 11.455e-6

        refl_mult_red = float(mtl.get("REFLECTANCE_MULT_BAND_3", "1"))
        refl_add_red = float(mtl.get("REFLECTANCE_ADD_BAND_3", "0"))
        refl_mult_nir = float(mtl.get("REFLECTANCE_MULT_BAND_4", "1"))
        refl_add_nir = float(mtl.get("REFLECTANCE_ADD_BAND_4", "0"))

    # Validate bands
    if not (red_path and nir_path and thr_path):
        raise arcpy.ExecuteError(f"Missing required bands in {scene_folder}")

    # Output naming
    scene_id = mtl.get("LANDSAT_SCENE_ID", os.path.basename(scene_folder))
    lst_path = out_lst_path

    ndvi_path = os.path.join(os.path.dirname(lst_path), f"NDVI_{scene_id}.tif") if save_ndvi else None
    emissivity_path = os.path.join(os.path.dirname(lst_path), f"EMIS_{scene_id}.tif") if save_emissivity else None
    bt_path = os.path.join(os.path.dirname(lst_path), f"BT_{scene_id}.tif") if save_bt else None

    arcpy.CheckOutExtension("Spatial")

    red = Raster(red_path)
    nir = Raster(nir_path)
    thr = Raster(thr_path)

    # ---- NDVI ----
    red_toa = red * refl_mult_red + refl_add_red
    nir_toa = nir * refl_mult_nir + refl_add_nir
    ndvi = (nir_toa - red_toa) / (nir_toa + red_toa)

    # ---- Pv ----
    pv = ((ndvi - 0.2) / (0.5 - 0.2)) ** 2
    pv = Con(pv < 0, 0, Con(pv > 1, 1, pv))

    # ---- Emissivity ----
    emissivity = 0.004 * pv + 0.986
    emissivity = Con(emissivity < 0.97, 0.97, Con(emissivity > 0.995, 0.995, emissivity))

    # ---- Radiance & BT ----
    radiance = thr * rad_mult + rad_add
    bt = k2 / Ln((k1 / radiance) + 1)

    # ---- LST ----
    rho = 1.4388e-2
    lst_k = bt / (1 + (lambda_thermal * bt / rho) * Ln(emissivity))
    lst_c = lst_k - 273.15

    lst_c.save(lst_path)

    if save_ndvi: ndvi.save(ndvi_path)
    if save_emissivity: emissivity.save(emissivity_path)
    if save_bt: bt.save(bt_path)

    arcpy.CheckInExtension("Spatial")

    return lst_path


# ----------------------------------------------------------
# Summary & Report
# ----------------------------------------------------------
def get_raster_stats(path):
    arcpy.management.CalculateStatistics(path)
    return {
        "min": float(arcpy.GetRasterProperties_management(path, "MINIMUM").getOutput(0)),
        "max": float(arcpy.GetRasterProperties_management(path, "MAXIMUM").getOutput(0)),
        "mean": float(arcpy.GetRasterProperties_management(path, "MEAN").getOutput(0)),
        "std": float(arcpy.GetRasterProperties_management(path, "STD").getOutput(0))
    }


def write_report(lst_path, stats, use_llm, api_key):
    folder = os.path.dirname(lst_path)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = os.path.join(folder, f"LST_report_{timestamp}.txt")

    with open(report_path, "w") as f:
        f.write(f"Land Surface Temperature Report\n")
        f.write(f"Raster: {os.path.basename(lst_path)}\n")
        f.write(f"Min: {stats['min']}\nMax: {stats['max']}\nMean: {stats['mean']}\nStd: {stats['std']}\n")

        if use_llm:
            response = call_llm(api_key, stats)
            f.write("\nLLM Interpretation:\n" + response if response else "\nNo LLM summary\n")

    return report_path


def call_llm(api_key, stats):
    if not api_key or not requests: return None

    prompt = f"""
    Interpret LST values:
    Min={stats['min']}, Max={stats['max']}, Mean={stats['mean']}, Std={stats['std']}
    """

    try:
        resp = requests.post(
            "https://api.openai.com/v1/chat/completions",
            headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
            json={"model": "gpt-4.1-mini", "messages": [{"role": "user", "content": prompt}]}
        )
        return resp.json()["choices"][0]["message"]["content"]

    except:
        return None


# ----------------------------------------------------------
# Main
# ----------------------------------------------------------
def main():
    arcpy.env.overwriteOutput = True

    scene_or_parent = arcpy.GetParameterAsText(0)
    output_folder = arcpy.GetParameterAsText(1)
    thermal_band_param = arcpy.GetParameter(2)
    save_ndvi = bool(arcpy.GetParameter(3))
    save_emissivity = bool(arcpy.GetParameter(4))
    save_bt = bool(arcpy.GetParameter(5))
    batch_mode = bool(arcpy.GetParameter(6))
    use_llm = bool(arcpy.GetParameter(7))
    api_key = get_api_key(arcpy.GetParameterAsText(8))

    thermal_band_number = int(thermal_band_param) if thermal_band_param else None

    if not os.path.isdir(scene_or_parent):
        log_error("Invalid folder.")
        return

    os.makedirs(output_folder, exist_ok=True)

    if not batch_mode:
        # ---- Single scene ----
        scene_id = os.path.basename(scene_or_parent)
        out_path = os.path.join(output_folder, f"LST_{scene_id}.tif")

        lst_path = compute_landsat_lst_for_scene(scene_or_parent, thermal_band_number, out_path,
                                                 save_ndvi, save_emissivity, save_bt)

        stats = get_raster_stats(lst_path)
        write_report(lst_path, stats, use_llm, api_key)

        arcpy.SetParameterAsText(1, lst_path)

    else:
        # ---- Batch Mode ----
        for name in os.listdir(scene_or_parent):
            subfolder = os.path.join(scene_or_parent, name)
            if not os.path.isdir(subfolder): continue

            try:
                scene_id = os.path.basename(subfolder)
                out_path = os.path.join(output_folder, f"LST_{scene_id}.tif")

                lst_path = compute_landsat_lst_for_scene(subfolder, thermal_band_number, out_path,
                                                         save_ndvi, save_emissivity, save_bt)

                stats = get_raster_stats(lst_path)
                write_report(lst_path, stats, use_llm, api_key)

                log_message(f"✔ Finished: {lst_path}")

            except Exception as e:
                log_warning(f"⚠ Skipped {subfolder} — {e}")


if __name__ == "__main__":
    main()

