#!/usr/bin/env python

import sys
import pandas as pd

intable = sys.argv[1]
outable=sys.argv[2]
#out_rute="GenePredAnnot/interproSCAN_outputs/parsed/"+outable

data= pd.read_table(intable,sep='\t')

g2 =data.groupby(['ID'])["Analysis"].apply(list).reset_index()
g3 =data.groupby(['ID'])["AnAccess"].apply(list).reset_index()
g4 =data.groupby(['ID'])["Description"].apply(list).reset_index()
g5 = data.groupby(['ID'])["InterPro_ann"].apply(list).reset_index()
g6 = data.groupby(['ID'])["Interpro_description"].apply(list).reset_index()
g7 = data.groupby(['ID'])["GO"].apply(list).reset_index()
g8 = data.groupby(['ID'])["Pathways"].apply(list).reset_index()
gN = data[["ID"]].drop_duplicates()

df2 = g2.merge(gN)
df2 = df2.merge(g3)
df2 = df2.merge(g4)
df2 = df2.merge(g5)
df2 = df2.merge(g6)
df2 = df2.merge(g7)
df2 = df2.merge(g8)



df2.to_csv(outable, index=None, sep='\t')

