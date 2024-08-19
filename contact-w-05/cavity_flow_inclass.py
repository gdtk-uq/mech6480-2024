"""
MECH6480 - WEEK 5 - Cavity Flow Example
We will develop this code in the Contact session on Wednesday.

Problem description:
 Re = 100
 rho = 1
 nu = 0.01
 Lx = Ly = 1
 
 u = Re*nu/Lx
   = 1e2*1e-2/1
   = 1
 
 Problem domain:
        lid: u = 1, v = 0
        +->->->->->->->->->->->->->-+
        |                           |
        |                           |
        |                           |
  wall  |                           |  wall
 u = 0  |                           |  u = 0
 v = 0  |                           |  v = 0
        |                           |
        |                           |
        +---------------------------+
        wall: u = 0, v = 0


 The staggered grid with ghost cells: 

    •   →   •   →   •   →   •   →   •
        |       |       |       |    
    ↑ - +---↑---+---↑---+---↑---+ - ↑
        :       |       |       :    
    •   →   •   →   •   →   •   →   •
        :       |       |       :    
    ↑ - + - ↑ - + - ↑ - + - ↑ - + - ↑
        :       |       |       :    
    •   →   •   →   •   →   •   →   •
        :       |       |       :    
    ↑ - + - ↑ - + - ↑ - + - ↑ - + - ↑
        :       |       |       :    
    •   →   •   →   •   →   •   →   •
        :       |       |       :    
    ↑ - 0---↑---+---↑---+---↑---+ - ↑
        |       |       |       |    
    •   →   •   →   •   →   •   →   •

 • Pressure stored at the cell centers
 → Horizontal velocity stored at the cell faces
 ↑ Vertical velocity stored at the cell faces
 0 Indicates origin of the grid
        
"""
