"""Rebuild the standard parked blades and digitize the measured Risoe curves."""
from pathlib import Path
import csv
import re
import struct
import numpy as np
import pdfplumber

root = Path(__file__).resolve().parent
with pdfplumber.open(root / 'nrel.pdf') as doc:
    plan = doc.pages[73].extract_text().splitlines()
    profile = doc.pages[74].extract_text().splitlines()
rows = []
for line in plan:
    if not re.match(r'^[0-5]\.\d+ ', line) or line.startswith('0.0 '):
        continue
    nums = re.findall(r'-?\d+(?:\.\d+)?', line)
    r = float(nums[0][:5])  # remove PDF footnote superscripts from radial distances
    if r == .508:
        rows.append([r,.218,0,.218,.5])
    else:
        chord, twist = map(float, nums[3:5])
        thick = .2095*chord if '%' in line else float(nums[5])
        rows.append([r,chord,twist,thick,float(nums[-1])/100])
plan = np.array(rows)
air = np.array([list(map(float,line.split())) for line in profile if re.fullmatch(r'[\d. -]+',line) and len(line.split())==4])
upper = np.vstack(([0,0],air[:,:2]))
lower = np.vstack(([0,0],air[:-1,2:]))
R = 5.029
def station(r):
    return [np.interp(r,plan[:,0],col) for col in plan.T[1:]]
def surface(r,x,side):
    c,twist,t,axis = station(r)
    theta = np.deg2rad(twist-station(R)[1])
    prof = upper if side==1 else lower
    b = np.interp(x,prof[:,0],prof[:,1])
    w = np.clip((r-.883)/(.374),0,1)
    h = ((1-w)*side*np.sqrt(np.maximum(0,x*(1-x))) + w*b/.2095)*t
    v = (x-axis)*c
    return np.array([np.sin(theta)*v+np.cos(theta)*h, np.cos(theta)*v-np.sin(theta)*h, np.full_like(x,r)]).T

# Root transition interpolates from the reported circular adapter to the S809 profile.
xs = np.unique(np.r_[upper[:,0],lower[:,0],0,1])
rings = []
for r in np.r_[plan[plan[:,0]<R,0],R]:
    rings.append(np.vstack((surface(r,xs,1),surface(r,xs[-2:0:-1],-1))))
triangles = []
for a,b in zip(rings,rings[1:]):
    for j in range(len(a)):
        k = (j+1)%len(a)
        triangles.extend([[a[j],a[k],b[k]],[a[j],b[k],b[j]]])
for ring,flip in [(rings[0],True),(rings[-1],False)]:
    for j in range(len(ring)):
        tri = [ring.mean(axis=0),ring[j],ring[(j+1)%len(ring)]]
        triangles.append(tri[::-1] if flip else tri)
triangles = np.array(triangles)
# Check orientation using the signed volume of the closed blade.
if np.einsum('ij,ij->i',triangles[:,0],np.cross(triangles[:,1],triangles[:,2])).sum()<0:
    triangles=triangles[:,::-1]
triangles=np.concatenate([triangles,triangles*np.array([1,-1,-1])])
with (root/'rotor.stl').open('wb') as f:
    f.write(bytes(80)+struct.pack('<I',len(triangles)))
    for t in triangles:
        n=np.cross(t[1]-t[0],t[2]-t[0])
        f.write(struct.pack('<12fH',*(n/np.linalg.norm(n)),*t.ravel(),0))

with pdfplumber.open(root/'risoe.pdf') as doc:
    page=doc.pages[20]
    loop=page.curves[79]['pts']  # Figure 15 r=47%, measured red trace with triangle markers
    cntrace=page.curves[0]['pts']  # Figure 14 measured sectional normal force
    assert len(loop)==23 and len(cntrace)==5
    # Determine grid lines from the vector plot, coordinates in PDF points from top.
    cp_zero=378.462
    cp_unit=34.600
    x0,x1=275.8116,403.5912
    points=[]
    for i,(x,y) in enumerate(loop[:-1]):
        xc=round((x-x0)/(x1-x0),3)
        side=1 if i<=12 else -1
        points.append((side,xc,-(cp_zero-y)/cp_unit))
    # Fig. 15 uses positive ordinate on the leeward side. Use the physical
    # pressure convention (p-p_inf)/q. The normal load gives an independent check.
    cn=(193.1552-cntrace[1][1])/138.4+.8
with (root/'probes.csv').open('w',newline='') as f:
    out=csv.writer(f)
    out.writerow(['side','xc','Cp_exp','X','Y','Z'])
    for side,xc,cp in points:
        out.writerow([side,xc,cp,*surface(.47*R,np.array([xc]),side)[0]])
    for side in (1,-1):
        for xc in np.linspace(0,1,51):
            out.writerow([side,xc,'nan',*surface(.47*R,np.array([xc]),side)[0]])
(root/'reference_cn.txt').write_text(f'{cn:.6f}\n')
np.savetxt(root/'planform.csv',plan,delimiter=',',header='r,chord,twist_deg,thickness,axis_fraction',comments='')
