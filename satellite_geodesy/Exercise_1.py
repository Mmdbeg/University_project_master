from dependencies_for_ex1 import *


# LOADING FILES PATH 
file_path_galileo = "galileo/GAL_only_nav.rnx"
file_path_GPS = "GPS/GPS_only_nav.rnx"
file_path_beidu = "beidu/BDS_only_nav.rnx"


#  GPS SATTELITE SYSTEM +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def main(path,rinex_content=""):
    try:

        with open(path, 'r', encoding='utf-8') as file:
            rinex_content = file.read()

        print(f"'{path}' loaded succesfully")


        if rinex_content:
            df_keplerian_multi = process_rinex_file_to_dataframe(rinex_content)


            #print("\n--- keplerian parammeters ---")
            #print(df_keplerian_multi)
        else:
            print("epmty file")

    except FileNotFoundError:
        print(f"file error '{path}' not found ")
    except Exception as e:
        print(f"an error accorded {e}")

    return df_keplerian_multi 


df_GPS = main(file_path_GPS)

df_GPS['Eccentric_Anomaly(new_raph)'] = None
df_GPS['Eccentric_Anomaly(expansion)'] = None

df_GPS['Eccentric_Anomaly(new_raph)'] = df_GPS.apply(
    lambda row: newton_raphson(row['M0 (rad)'], row['e'], row['M0 (rad)']), 
    axis=1
)
df_GPS['Eccentric_Anomaly(expansion)'] = df_GPS.apply(
    lambda row: series_eccentric_anomaly(row['M0 (rad)'], row['e']), 
    axis=1
)




#  BEIDU SATTELITE SYSTEM +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def main(path,rinex_content=""):
    try:

        with open(path, 'r', encoding='utf-8') as file:
            rinex_content = file.read()

        print(f"'{path}' loaded succesfully")


        if rinex_content:
            df_keplerian_multi = process_rinex_file_to_dataframe(rinex_content)


            #print("\n--- keplerian parammeters ---")
            #print(df_keplerian_multi)
        else:
            print("epmty file")

    except FileNotFoundError:
        print(f"file error '{path}' not found ")
    except Exception as e:
        print(f"an error accorded {e}")

    return df_keplerian_multi 


df_beidu = main(file_path_beidu)

df_beidu['Eccentric_Anomaly(new_raph)'] = None
df_beidu['Eccentric_Anomaly(expansion)'] = None

df_beidu['Eccentric_Anomaly(new_raph)'] = df_beidu.apply(
    lambda row: newton_raphson(row['M0 (rad)'], row['e'], row['M0 (rad)']), 
    axis=1
)
df_beidu['Eccentric_Anomaly(expansion)'] = df_beidu.apply(
    lambda row: series_eccentric_anomaly(row['M0 (rad)'], row['e']), 
    axis=1
)





#   GALILEO SATTELITE SYSTEM +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def main(path,rinex_content=""):
    try:

        with open(path, 'r', encoding='utf-8') as file:
            rinex_content = file.read()

        print(f"'{path}' loaded succesfully")


        if rinex_content:
            df_keplerian_multi = process_rinex_file_to_dataframe(rinex_content)


            #print("\n--- keplerian parammeters ---")
            #print(df_keplerian_multi)
        else:
            print("epmty file")

    except FileNotFoundError:
        print(f"file error '{path}' not found ")
    except Exception as e:
        print(f"an error accorded {e}")

    return df_keplerian_multi 


df_galileo = main(file_path_galileo)

df_galileo['Eccentric_Anomaly(new_raph)'] =None
df_galileo['Eccentric_Anomaly(expansion)'] =None


df_galileo['Eccentric_Anomaly(new_raph)'] = df_galileo.apply(
    lambda row: newton_raphson(row['M0 (rad)'], row['e'], row['M0 (rad)']), 
    axis=1
)
df_galileo['Eccentric_Anomaly(expansion)'] = df_galileo.apply(
    lambda row: series_eccentric_anomaly(row['M0 (rad)'], row['e']), 
    axis=1
)




print(df_GPS)
print(df_beidu)
print(df_galileo)

