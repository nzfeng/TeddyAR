
 using UnityEngine;
 
 public static class Vector3Extensions
 {
     public static Vector2 XY(this Vector3 aVector)
     {
         return new Vector2(aVector.x,aVector.y);
     }
     public static Vector2 XZ(this Vector3 aVector)
     {
         return new Vector2(aVector.x,aVector.z);
     }
     public static Vector2 YZ(this Vector3 aVector)
     {
         return new Vector2(aVector.y,aVector.z);
     }
     public static Vector2 YX(this Vector3 aVector)
     {
         return new Vector2(aVector.y,aVector.x);
     }
     public static Vector2 ZX(this Vector3 aVector)
     {
         return new Vector2(aVector.z,aVector.x);
     }
     public static Vector2 ZY(this Vector3 aVector)
     {
         return new Vector2(aVector.z,aVector.y);
     }
 }