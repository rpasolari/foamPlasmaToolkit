# Tutorial: Implicit Boundary Coupling (fluid + dielectric)

**Purpose**  
This tutorial demonstrates how **implicit boundary coupling** works between two regions (a *fluid*
and a *dielectric*) — using the `coupledElectricPotential` boundary condition.  
In this mode, the interface is treated **implicitly**, meaning a **single global linear system** is 
assembled and solved **monolithically across both regions**.

## 1. How to Run

You can run the tutorial in **serial** or **parallel**:

```bash
# Serial execution
./Allrun-serial

# Parallel execution
./Allrun-parallel
```

>[!IMPORTANT]
>For parallel runs, it is essential that the owner and neighbour cells of every face on the 
interface patch lie within the same processor. Otherwise, the implicit coupling will fail because 
the global matrix cannot properly connect both sides of the interface. 
Check your `system/decomposeParDict`.

