from functions import *

data_path = "data.xlsx"

df = pd.read_excel(data_path)
pd.set_option('display.float_format', '{:.10f}'.format)

# Create empty columns for results
df["Longitude"] = None
df["Latitude"] = None

# you can use py-proj instead  of utm to latlon function (that developed by programmer) un-comment this part VVVV
'''
# for i, row in df.iterrows():
#     zone = int(row["Zone"])        # read zone from your data
#     x = row["E"]
#     y = row["N"]

#     # Iran is in the northern hemisphere
#     transformer = pj.Transformer.from_crs(
#         f"+proj=utm +zone={zone} +ellps=WGS84 +north",
#         "EPSG:4326",  # WGS84 geographic
#         always_xy=True
#     )

#     lon, lat = transformer.transform(x, y)
#     df.at[i, "Longitude"] = lon
#     df.at[i, "Latitude"] = lat

# ^^^^^^^^^^^^^^^
'''

# Assume all are in northern hemisphere (Iran)
df["Latitude"], df["Longitude"] = zip(*df.apply(lambda r: utm_to_latlon(r["E"], r["N"], int(r["Zone"]), northern=True), axis=1))


#print(df[["E", "N", "Zone", "Longitude" , "Latitude"]])

# +-+-+-+-+-+- END OF PART 2 :  UTM (X , Y) > WGS-84 ( LAT , LONG )  +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

geoid = GeoidPGM('GeographicLib/geoids/egm96-5.pgm')  


for i in range(len(df)):

    lat = index(df,i,"Latitude")
    lon = index(df,i,"Longitude")
    H   = index(df,i,"H")   # Google Earth orthometric height H

    # 1. Geoid undulation N
    N = geoid.height(lat, lon)
    df.loc[i,'geoidal_undulation'] = N

    # 2. Ellipsoidal height
    h = H + N
    df.loc[i,'h'] = h

    # 3. Convert (lat, lon, h) → geocentric XYZ
    l1 = LatLon(lat, lon, height=h)
    cartesian_obj = l1.toCartesian()

    df.loc[i,"X_wgs84"] = cartesian_obj.x
    df.loc[i,"Y_wgs84"] = cartesian_obj.y
    df.loc[i,"Z_wgs84"] = cartesian_obj.z

    #print(index(df,i,"h") ,"-->" , index(df,i,"geoidal_undulation"))



# Extracting geocentre cartesian coordinates in WGS84
xyz = df[df.columns[-3:]]
xyz = xyz.copy()  # now xyz is independent

for i in range(len(xyz)):
    x, y, z = wgs84_to_itrf2000(
        xyz.loc[i, "X_wgs84"],
        xyz.loc[i, "Y_wgs84"],
        xyz.loc[i, "Z_wgs84"]
    )
    xyz.loc[i, "X_ITRF2000"] = x
    xyz.loc[i, "Y_ITRF2000"] = y
    xyz.loc[i, "Z_ITRF2000"] = z
    
# WGS84(LAT,LON,h >> x,y,z) to ITRF2000 (X,Y,X) END OF PART 3 +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-    

# for i in range(len(xyz)):    
#     print(np.sqrt( (xyz.loc[i, "X_ITRF2000"]-xyz.loc[i, "X_wgs84"])**2 + (xyz.loc[i, "Y_ITRF2000"]-xyz.loc[i, "Y_wgs84"])**2 + (xyz.loc[i, "Z_ITRF2000"]-xyz.loc[i, "Z_wgs84"])**2 ))

# 1. converting coordinates XYZ (ITRF2000)  --->  XYZ (ITRF2005)   //// IN refrence epoch (2000) ////

# displacement vector 
r1 =[]

#  converting ITRF2000 --> ITRF2005 (EPOCH 2000)
for i in range(len(xyz)):

    x1, y1, z1 = itrf2000_to_itrf2005(
        xyz.loc[i, "X_ITRF2000"], 
        xyz.loc[i, "Y_ITRF2000"],
        xyz.loc[i, "Z_ITRF2000"] 
    )
    xyz.loc[i, "X_ITRF2005"] = x1
    xyz.loc[i, "Y_ITRF2005"] = y1
    xyz.loc[i, "Z_ITRF2005"] = z1

    r1.append(np.sqrt( (xyz.loc[i, "X_ITRF2000"]-x1)**2 + (xyz.loc[i, "Y_ITRF2000"]-y1)**2 + (xyz.loc[i, "Z_ITRF2000"]-z1)**2 ))

it2000to2005_displacement = np.asarray(r1) 

# 2. converting coordinates XYZ (ITRF2005)  --->  XYZ (ITRF2008)   //// IN refrence epoch (2000) ////
r2 = [ ]
for i in range(len(xyz)):

    x2, y2, z2 = itrf2005_to_itrf2008(
        xyz.loc[i, "X_ITRF2005"], 
        xyz.loc[i, "Y_ITRF2005"],
        xyz.loc[i, "Z_ITRF2005"] 
    )
    xyz.loc[i, "X_ITRF2008"] = x2
    xyz.loc[i, "Y_ITRF2008"] = y2
    xyz.loc[i, "Z_ITRF2008"] = z2

    r2.append(np.sqrt( (xyz.loc[i, "X_ITRF2005"]-x2)**2 + (xyz.loc[i, "Y_ITRF2005"]-y2)**2 + (xyz.loc[i, "Z_ITRF2005"]-z2)**2 ))

it2005to2008_displacement = np.asarray(r2) 


# 3. converting coordinates XYZ (ITRF2008)  --->  XYZ (ITRF2014)   //// IN refrence epoch (2010) ////
r3 = [ ]
for i in range(len(xyz)):

    x3, y3, z3 = itrf2008_to_itrf2014(
        xyz.loc[i, "X_ITRF2008"], 
        xyz.loc[i, "Y_ITRF2008"],
        xyz.loc[i, "Z_ITRF2008"] 
    )
    xyz.loc[i, "X_ITRF2014"] = x3
    xyz.loc[i, "Y_ITRF2014"] = y3
    xyz.loc[i, "Z_ITRF2014"] = z3

    r3.append(np.sqrt( (xyz.loc[i, "X_ITRF2008"]-x3)**2 + (xyz.loc[i, "Y_ITRF2008"]-y3)**2 + (xyz.loc[i, "Z_ITRF2008"]-z3)**2 ))

it2008to2014_displacement = np.asarray(r3) 

# 4. converting coordinates XYZ (ITRF2014)  --->  XYZ (ITRF2020)   //// IN refrence epoch (2015) ////
r4 = [ ]
for i in range(len(xyz)):

    x4, y4, z4 = itrf2014_to_itrf2020(
        xyz.loc[i, "X_ITRF2014"], 
        xyz.loc[i, "Y_ITRF2014"],
        xyz.loc[i, "Z_ITRF2014"] 
    )
    xyz.loc[i, "X_ITRF2020"] = x4
    xyz.loc[i, "Y_ITRF2020"] = y4
    xyz.loc[i, "Z_ITRF2020"] = z4

    r4.append(np.sqrt( (xyz.loc[i, "X_ITRF2014"]-x4)**2 + (xyz.loc[i, "Y_ITRF2014"]-y4)**2 + (xyz.loc[i, "Z_ITRF2014"]-z4)**2 ))

it2014to2020_displacement = np.asarray(r4) 


print("it2000to2005_displacement :",it2000to2005_displacement,"\n")
print("it2005to2008_displacement :",it2005to2008_displacement,"\n")
print("it2008to2014_displacement :",it2008to2014_displacement,"\n")
print("it2014to2020_displacement :",it2014to2020_displacement,"\n")
print("****************************************************************")
# #------------------------------------------------------------------------------------------------------------------------------

# # PART 6 ---- ITRF2000 --> ITRF2005 (EPOCH 2022)  ////// ITRF2000 --> ITRF2022  (EPOCH 2022)  

# a subset of 10 points +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
xyz_10 = xyz.head(10)
xyz_10 = xyz_10.copy()  # now xyz_10 is independent

# converting 10 points from itrf2000 to itrf2020 (in 2022 epoch)

r_subset_2000_2005 = [ ]
r_subset_2000_2020 = [ ]


for i in range(len(xyz_10)):

    xs2, ys2, zs2 = itrf2000_to_itrf2020_in2022(
        xyz_10.loc[i, "X_ITRF2000"], 
        xyz_10.loc[i, "Y_ITRF2000"],
        xyz_10.loc[i, "Z_ITRF2000"] 
    )
    xyz_10.loc[i, "X_ITRF2020_2022"] = xs2
    xyz_10.loc[i, "Y_ITRF2020_2022"] = ys2
    xyz_10.loc[i, "Z_ITRF2020_2022"] = zs2

    r_subset_2000_2020.append(np.sqrt( (xyz_10.loc[i, "X_ITRF2000"]-xs2)**2 + (xyz_10.loc[i, "Y_ITRF2000"]-ys2)**2 + (xyz_10.loc[i, "Z_ITRF2000"]-zs2)**2 ))

# converting 10 points from itrf2000 to itrf2005 (in 2022 epoch)

for i in range(len(xyz_10)):

    xs1, ys1, zs1 = itrf2000_to_itrf2005_in2022(
        xyz_10.loc[i, "X_ITRF2000"], 
        xyz_10.loc[i, "Y_ITRF2000"],
        xyz_10.loc[i, "Z_ITRF2000"] 
    )
    xyz_10.loc[i, "X_ITRF2005_2022"] = xs1
    xyz_10.loc[i, "Y_ITRF2005_2022"] = ys1
    xyz_10.loc[i, "Z_ITRF2005_2022"] = zs1

    r_subset_2000_2005.append(np.sqrt( (xyz_10.loc[i, "X_ITRF2000"]-xs1)**2 + (xyz_10.loc[i, "Y_ITRF2000"]-ys1)**2 + (xyz_10.loc[i, "Z_ITRF2000"]-zs1)**2 ))





it2000to2020_in2022__displacement = np.asarray(r_subset_2000_2020) 
it2000to2005_in2022__displacement = np.asarray(r_subset_2000_2005) 
print("it2000to2020_in2022__displacement :",it2000to2020_in2022__displacement,"\n")
print("it2000to2005_in2022__displacement : ",it2000to2005_in2022__displacement,"\n")





# Stack all vectors together to create a 2D array
data = np.array([it2000to2005_displacement, 
                 it2005to2008_displacement, 
                 it2008to2014_displacement, 
                 it2014to2020_displacement]).T  # Shape will be (20, 4)

# Labels for each time period
time_labels = ['2000-2005', '2005-2008', '2008-2014', '2014-2020']

# Position of the bars on the x-axis
x = np.arange(20)  # 20 points

# Bar width
bar_width = 0.2

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8))

# Create bars for each time period
for i, (time_label, displacement) in enumerate(zip(time_labels, data.T)):
    ax.bar(x + i * bar_width, displacement, bar_width, label=time_label)

# Customize the plot
ax.set_xlabel('Points')
ax.set_ylabel('Deformation')
ax.set_title('Deformation Comparison Across Time Periods')
ax.set_xticks(x + 1.5 * bar_width)
ax.set_xticklabels([f'p {i+1}' for i in range(20)])
ax.legend(title='Time Period')

# Show the plot
plt.tight_layout()
plt.show()