"""Оверлей Tier 0 (силы) на оригинальную мастер-фигуру Part 3.
Оригинал сохраняется пиксель-в-пиксель; добавляется ряд ранга 1 + сноска на расширении холста."""
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

fig=plt.figure(figsize=(W/100,(H+EXT)/100),dpi=100)
ax=fig.add_axes([0,0,1,1]); ax.imshow(canvas); ax.set_xlim(0,W); ax.set_ylim(H+EXT,0); ax.axis('off')

gold='#8a5a00'
y1=1194; x0,x1=470,1956
# линия ряда
ax.plot([x0,x1],[y1,y1],color=gold,lw=5,solid_capstyle='round',zorder=3)
# стартовый кластер: три отделения (слипшиеся на этой шкале)
for dx,ms in [(0,17),(26,15),(52,13)]:
    ax.plot([x0+dx],[y1],marker='v',ms=ms,color=gold,mec='black',mew=1.2,zorder=5)
# итоговый маркер в стиле колец + число
ax.plot([x1],[y1],marker='o',ms=34,color=gold,mec='black',mew=2.5,zorder=5)
ax.annotate('4',(x1+44,y1),fontsize=23,color=gold,fontweight='bold',va='center')
# подпись ряда (в стиле остальных)
ax.annotate('Forces (fundamental interactions)',(x0+85,y1-26),fontsize=21,color=gold,fontweight='bold')
# выноска-бокс в стиле рисунка
txt=('Forces differentiate ~13.8 Gyr ago\n(within the first 10⁻¹² s):\n'
     'gravity → strong (interference troughs:\nfirst composites) → EM / weak split\n(gradient / template relaxation)')
ax.annotate(txt,(x0+26,y1),xytext=(1480,1030),fontsize=15.5,color=gold,ha='center',va='top' if False else 'top',
            linespacing=1.35,
            bbox=dict(boxstyle='round,pad=0.45',fc='white',ec=gold,lw=1.6),
            zorder=6)
# сноска на расширении
ax.annotate('Overlay (this work): Tier 0 — differentiation of forces. Schematic: the entire sequence sits at ≈13.8 Gyr before present on this axis; the reading of each step follows Sections 3–4.',
            (262,H+44),fontsize=14.5,color=gold,style='italic')
plt.savefig('figure1_with_forces_tier0.png',dpi=100)
print("saved figure1_with_forces_tier0.png")
