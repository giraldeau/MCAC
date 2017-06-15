
#include "test-pov.inc"  
#include "colors.inc" 
#include "textures.inc"
#include "stones.inc"   

#declare L =319.95837; 
 
global_settings { ambient_light rgb<1, 1, 1>
  max_trace_level 10  }

camera { 
        orthographic 
   location  <L/1.5, L/1.5, -1.0*L>   
   //location  <0, 0, -2.0*L>  
   look_at   <0, 0, 0>
        }       
    
 plane {             
  z, L    
   texture {  pigment { White }   finish { ambient 0.4 }
   }  
   no_shadow                   
       }   
       
       
 difference {
   box { <-L/2, -L/2, -L/2>, <L/2, L/2, L/2> } 
   box { <-L/2*0.9, -L/2*0.9, -L/2*1.1>, <L/2*0.9, L/2*0.9, L/2*1.1> } 
   box { <-L/2*1.1, -L/2*0.9, -L/2*0.9>, <L/2*1.1, L/2*0.9, L/2*0.9> }
   box { <-L/2*0.9, -L/2*1.1, -L/2*0.9>, <L/2*0.9, L/2*1.1, L/2*0.9> }
   pigment {color rgbt <1.0, 0., 0.,0.>}    
   no_shadow    
   }
    
       
       
 
object { agreg    
 //translate <-L/2,-L/2,-L/2>
 rotate <0, 0, 0>
 pigment {color rgbt <0.1, 0.1, 0.1,0.>}     
 normal { granite 0.2 }
 normal { bumps 1 scale 3 }   
 no_shadow                

//interior{ ior 1.8 }
 //finish { reflection {0.0} specular 0.5 roughness 0.9 ambient 0.0}  
 }

light_source { <0, 0, -50000> color red 1 green 1 blue 1}
