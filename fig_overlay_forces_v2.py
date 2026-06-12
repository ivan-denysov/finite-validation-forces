"""Оверлей v2: правая панель убрана; вместо неё — силы, каждая со стрелкой-нитью
к старту кольца, которое она открывает. Ряд Tier 0 (когда отделились) сохранён."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

src=np.array(Image.open('/mnt/user-data/uploads/figure1_four_ring_continuum.png'))
H,W,_=src.shape
EXT=70
canvas=np.full((H+EXT,W,3),255,dtype=np.uint8)
canvas[:H]=src[:,:,:3]
canvas[140:1420,2210:W]=255   # убрать правую панель

fig=plt.figure(figsize=(W/100,(H+EXT)/100),dpi=100)
ax=fig.add_axes([0,0,1,1]); ax.imshow(canvas); ax.set_xlim(0,W); ax.set_ylim(H+EXT,0); ax.axis('off')
gold='#8a5a00'

# ---- Tier 0 ряд (как в v1) ----
y1=1194; x0,x1=470,1956
ax.plot([x0,x1],[y1,y1],color=gold,lw=5,solid_capstyle='round',zorder=3)
for dx,ms in [(0,17),(26,15),(52,13)]:
    ax.plot([x0+dx],[y1],marker='v',ms=ms,color=gold,mec='black',mew=1.2,zorder=5)
ax.plot([x1],[y1],marker='o',ms=34,color=gold,mec='black',mew=2.5,zorder=5)
ax.annotate('4',(x1+44,y1),fontsize=23,color=gold,fontweight='bold',va='center')
ax.annotate('Forces (fundamental interactions)',(x0+85,y1-26),fontsize=21,color=gold,fontweight='bold')

# ---- Правая колонка: силы открывают кольца ----
PX=2245
ax.annotate('Forces open the rings',(PX,235),fontsize=26,color='black',fontweight='bold')
ax.annotate('(overlay, this work): each new ring arrives\ntogether with its own type of binding',
            (PX,275),fontsize=16,color='#555555',style='italic',va='top',linespacing=1.3)
ax.annotate('black = established literature   ·   gold = this work (overlay)',
            (PX,345),fontsize=13.5,color='#777777',style='italic')

entries=[
 ('1','GRAVITY','the delay gradient collapses gas\ninto the first stars', 905, (268,1118)),
 ('2','STRONG','interference troughs bind nucleons\n(candidate, this work):\nfirst composite templates → elements', 730, (300,1010)),
 ('3','ELECTROMAGNETIC','chemical bonding: charges and\nphotons → molecules and minerals', 555, (360,935)),
 ('4','WEAK','template relaxation sets slow stellar\nburning (pp chain): billions of years\nof steady sunlight → time for life', 380, (380,745)),
]
for num,name,desc,ty,(tx,tyy) in entries:
    ax.plot([PX+18],[ty+6],marker='o',ms=21,color=gold,mec='black',mew=1.5,zorder=7)
    ax.annotate(num,(PX+18,ty+7),fontsize=17,color='white',fontweight='bold',ha='center',va='center',zorder=8)
    ax.annotate(name,(PX+48,ty+15),fontsize=20,color=gold,fontweight='bold',va='center')
    ax.annotate(desc,(PX+48,ty+38),fontsize=15.5,color='#333333',va='top',linespacing=1.3)
    ax.plot([tx],[tyy],marker='o',ms=24,color=gold,mec='black',mew=1.5,zorder=7,alpha=0.95)
    ax.annotate(num,(tx,tyy+1),fontsize=18,color='white',fontweight='bold',ha='center',va='center',zorder=8)
ax.annotate('numbers = order of differentiation;\nposition = the ring each force opens',
            (PX,1100),fontsize=15,color='#555555',style='italic',va='top',linespacing=1.35)
ax.annotate('Overlay (this work): Tier 0 — differentiation of forces (bottom row; the entire sequence sits at ≈13.8 Gyr on this axis). Right column: the ring each force opens (Sections 3–4).',
            (262,H+44),fontsize=14.5,color=gold,style='italic')
plt.savefig('figure1_with_forces_tier0.png',dpi=100)
print("saved")
