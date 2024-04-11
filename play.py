
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')

import numpy as np
import matplotlib.pyplot as plt



def names_to_df(names, df, cols=['Name', "D", "P", "H", "viscosity"]):
    idx_list = []
    for n in names:
        #print(n)
        dfbool = df['Name'].str.fullmatch(n)
        idx = int(float(np.argwhere(np.array(dfbool))))
        idx_list.append(idx)
        selection = df.iloc[idx_list, :][cols]

    return selection


def get_distance(df,name,names):


    testsolvent = names_to_df([name], df)
    Dt,Pt,Ht = testsolvent[['D','P','H']].to_numpy().ravel().tolist()
    dists,directs = [],[]
    for n in names:
        listsolvent = names_to_df([n],df)
        Dl, Pl, Hl = listsolvent[['D', 'P', 'H']].to_numpy().ravel().tolist()
        # Calculate the Euclidean distance).
        distance = (
            4*(Dl-Dt)**2+
            (Pl-Pt)**2+
            (Hl-Ht)**2)**0.5

        # calculate direction vector tuple
        resid = (Dl-Dt),(Pl-Pt),(Hl-Ht)
        resid = np.array(resid) / np.linalg.norm(resid)

        theta = np.degrees(np.arccos(resid[2])) # polar angle 0 above 90 horizontal 180 downwards
        phi = np.degrees(np.arctan2(resid[1],resid[0])) # azimutal angle, 0 along D, 90 along P, 180 -D
        #print(phi,theta)
        dists.append(distance)
        directs.append([phi,theta])

    return dists,directs



df = pd.read_csv('db.csv')



names = ["1-Butanol",
         #"Dimethylformamide",
         #"Dimethyl sulfoxide",
         #"Water",
         "2-Propanol",
         "Toluene",
         "Anisole",
         "Chlorobenzene",
         "Ethanol",
         "m-Xylene",
         "Mesitylene",
         "Diethyl ether",
         "Chloroform",
         "Ethyl acetate",
         "Butyl acetate",
         "o-dichlorobenzene",
         ]


names_solvs = ["Dimethylformamide",
         "Dimethyl sulfoxide",
               "Water"
         ]

names_type1 = [
        "1-Butanol",
         "2-Propanol",
         "Ethanol",
         ]


names_type2 = [
         "Toluene",
         "Chlorobenzene",
         "Chloroform",
         "Ethyl acetate",
         "Butyl acetate",
         "o-dichlorobenzene",
         ]


names_type3 = [
         "m-Xylene",
         "Mesitylene",
         "Diethyl ether",
         ]

names_miscsolv = ['Acetonitrile',
              'Ethylene glycol monoethyl ether', # 2Me
              ]

names_misccoord = [
              'Methyl-2-pyrrolidone', # NMP
              'Gamma butyrolactone', # GBL
              ]

names_merge = names_type1+names_type2+names_type3+names_solvs+names_miscsolv+names_misccoord

# todo check doubles, could use dicts keys

# todo add  tbp

dists,directs = get_distance(df,'Dimethylformamide',names_type1)


dists,directs = get_distance(df,'Dimethylformamide',names_type2)
dists,directs = get_distance(df,'Dimethylformamide',names_type3)



fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.set_aspect('equal')
ax.set_box_aspect((2, 1, 1))

fig.suptitle('3D HSP Plot')

typ1 = names_to_df(names_type1,df)
typ2 = names_to_df(names_type2,df)
typ3 = names_to_df(names_type3,df)
solvs = names_to_df(names_solvs,df)
miscsolv = names_to_df(names_miscsolv,df)
misccoord = names_to_df(names_misccoord,df)

coord_merge = names_to_df(names_merge,df)

s=100
ax.scatter(typ1['D'],typ1['P'],typ1['H'],color='blue',label='Type1',s=s)
ax.scatter(typ2['D'],typ2['P'],typ2['H'],color='green',label='Type2',s=s)
ax.scatter(typ3['D'],typ3['P'],typ3['H'],color='yellow',label='Type3',s=s)
ax.scatter(solvs['D'],solvs['P'],solvs['H'],color='black',label='DMF/DMSO',s=s)
ax.scatter(miscsolv['D'],miscsolv['P'],miscsolv['H'],color='brown',label='solvs',s=s)
ax.scatter(misccoord['D'],misccoord['P'],misccoord['H'],color='purple',label='coord',s=s)

for idx,name in enumerate(names_merge):
    print(coord_merge['D'].iloc[idx])
    ax.text(coord_merge['D'].iloc[idx],coord_merge['P'].iloc[idx],coord_merge['H'].iloc[idx], name, size=8, zorder=99,color='k')

ax.set_xlabel('D')
ax.set_ylabel('P')
ax.set_zlabel('H')
ax.zaxis.labelpad=-2

delta=20
D_min,P_min,H_min=10,0,0

ax.set_xlim(D_min,D_min+delta)
ax.set_ylim(P_min,P_min+delta)
ax.set_zlim(H_min,H_min+delta)

plt.legend()

plt.show()






#instan.grid[instan.grid['synonyms'].str.contains('propan').astype(bool)]



#instan.grid[instan.grid['synonyms'].str.contains('butanol').astype(bool)][['Name',"D","P","H","viscosity"]]




