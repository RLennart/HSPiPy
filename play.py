import hsp as h

instan = h.HSP()

instan.read('db.csv')

instan.grid

#instan.grid[['Solvent','D','P','H','Score']]

instan.get(3)

