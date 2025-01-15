[AC Analysis]
{
   Npanes: 3
   {
      traces: 1 {2,0,"V(vout)/(V(vinp)-V(vinn))"}
      X: ('G',4,1000,0,1e+11)
      Y[0]: (' ',0,9.99999999999999e-05,10,316.227766016838)
      Y[1]: (' ',0,-330,30,60)
      Log: 1 2 0
      GridStyle: 1
      PltMag: 1
      PltPhi: 1 0
   },
   {
      traces: 1 {3,0,"V(vinter)/(V(vinn)-V(vinp))"}
      X: ('G',4,1000,0,1e+11)
      Y[0]: (' ',0,0.0630957344480193,4,15.8489319246111)
      Y[1]: (' ',0,-200,20,40)
      Log: 1 2 0
      GridStyle: 1
      PltMag: 1
      PltPhi: 1 0
   },
   {
      traces: 1 {4,0,"V(vout)/V(vinter)"}
      X: ('G',4,1000,0,1e+11)
      Y[0]: (' ',0,0.000630957344480193,8,100)
      Y[1]: (' ',0,50,10,180)
      Log: 1 2 0
      GridStyle: 1
      PltMag: 1
      PltPhi: 1 0
   }
}
[Transient Analysis]
{
   Npanes: 3
   Active Pane: 1
   {
      traces: 1 {524290,0,"V(vout)/(V(vinp)-V(vinn))"}
      X: (' ',1,0,0.1,1.024)
      Y[0]: ('_',0,-1e+307,1e+307,1e+308)
      Y[1]: (' ',0,1e+308,30,-1e+308)
      Units: "" ('_',0,0,1,-1e+307,1e+307,1e+308)
      Log: 0 0 0
      GridStyle: 1
      PltMag: 1
      PltPhi: 1 0
   },
   {
      traces: 1 {268959749,0,"V(vout)"}
      X: (' ',1,0,0.1,1.024)
      Y[0]: ('m',0,0,0.09,0.99)
      Y[1]: ('K',2,1e+308,60,-1e+308)
      Volts: ('m',0,0,1,0,0.09,0.99)
      Log: 0 0 0
      GridStyle: 1
      PltMag: 1
      PltPhi: 1 0
   },
   {
      traces: 1 {524292,0,"V(vout)/V(vinter)"}
      X: (' ',1,0,0.1,1.024)
      Y[0]: (' ',1,0,0.4,4.8)
      Y[1]: (' ',0,1e+308,20,-1e+308)
      Units: "" (' ',0,0,1,0,0.4,4.8)
      Log: 0 0 0
      GridStyle: 1
      PltMag: 1
      PltPhi: 1 0
   }
}
