#include "mol.h"

/*
  This function should copy the values pointed to by element, x, y, and z into the atom stored at
  atom. You may assume that sufficient memory has been allocated at all pointer addresses.
  Note that using pointers for the function “inputs”, x, y, and z, is done here to match the
  function arguments of atomget.
*/
void atomset(atom *atom, char element[3], double *x, double *y, double *z ){
  atom->x = *x;
  atom->y = *y;
  atom->z = *z;
  strcpy(atom->element, element);
}

/*
  This function should copy the values in the atom stored at atom to the locations pointed to by
  element, x, y, and z. You may assume that sufficient memory has been allocated at all pointer
  addresses. Note that using pointers for the function “input”, atom, is done here to match the
  function arguments of atomset.
*/
void atomget(atom *atom, char element[3], double *x, double *y, double *z){
  *x = atom->x;
  *y = atom->y;
  *z = atom->z;
  strcpy(element, atom->element);
}

/*
    This function should copy the values pointed to by a1, a2, atoms, and epairs into the 
  corresponding structure attributes in bond. In addition, you should call the compute_coords
  function (see below) on the bond.
*/
void bondset( bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs ){
  bond->a1 = *a1;
  bond->a2 = *a2;
  bond->atoms = *atoms;
  bond->epairs = *epairs;
  compute_coords(bond);
}

/*
    This function should copy the structure attributes in bond to their corresponding arguments:
  a1, a2 and epairs. You may assume that sufficient memory has been allocated at all pointer
  addresses. Note you are not copying atom structures, only the addresses of the atom
  structures.
*/
void bondget( bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs ){
  *a1 = bond->a1;
  *a2 = bond->a2;
  *epairs = bond->epairs;
  atoms = &bond->atoms;
}

/*
    This function should compute the z, x1, y1, x2, y2, len, dx, and dy values of the bond and set 
  them in the appropriate structure member variables.
*/
void compute_coords( bond *bond ){
  bond->x1 = bond->atoms[bond->a1].x;
  bond->x2 = bond->atoms[bond->a2].x;
  bond->y1 = bond->atoms[bond->a1].y;
  bond->y2 = bond->atoms[bond->a2].y;
  bond->z = (bond->atoms[bond->a1].z + bond->atoms[bond->a2].z) / 2;
  bond->len = pow(pow((bond->atoms[bond->a2].x - bond->atoms[bond->a1].x), 2) + pow((bond->atoms[bond->a2].y - bond->atoms[bond->a2].y), 2) + pow((bond->atoms[bond->a1].z - bond->atoms[bond->a2].z), 2), (1/2));
  bond->dx = (bond->atoms[bond->a2].x - bond->atoms[bond->a1].x) / bond->len;
  bond->dy = (bond->atoms[bond->a2].y - bond->atoms[bond->a1].y) / bond->len;
}

/*
  This function should return the address of a malloced area of memory, large enough to hold a
  molecule. The value of atom_max should be copied into the structure; the value of atom_no in
  the structure should be set to zero; and, the arrays atoms and atom_ptrs should be malloced
  to have enough memory to hold atom_max atoms and pointers (respectively). The value of
  bond_max should be copied into the structure; the value of bond_no in the structure should be
  set to zero; and, the arrays bonds and bond_ptrs should be malloced to have enough memory
  to hold bond_max bonds and pointers (respectively).
*/
molecule *molmalloc(unsigned short atom_max, unsigned short bond_max){
  molecule *new = malloc(sizeof(molecule));
  new->atom_max = atom_max;
  new->atom_no = 0;
  new->atoms = malloc(sizeof(atom) * atom_max);
  new->atom_ptrs = malloc(sizeof(atom*) * atom_max);

  new->bond_max = bond_max;
  new->bond_no = 0;
  new->bonds = malloc(sizeof(bond) * bond_max);
  new->bond_ptrs = malloc(sizeof(bond*) * bond_max);
  return new;
}

/*
  Same as molmalloc but uses realloc to resize.
*/
void atom_realloc(molecule *mol){
  mol->atoms = realloc(mol->atoms, sizeof(atom) * mol->atom_max);
  mol->atom_ptrs = realloc(mol->atom_ptrs, sizeof(atom*) * mol->atom_max);
}
//Reallocing for bonds
void bond_realloc(molecule *mol){
  mol->bonds = realloc(mol->bonds, sizeof(bond) * mol->bond_max);
  mol->bond_ptrs = realloc(mol->bond_ptrs, sizeof(bond*) * mol->bond_max);
}

/*
  This function should return the address of a malloced area of memory, large enough to hold a
  molecule. Additionally, the values of atom_max, atom_no, bond_max, bond_no should be
  copied from src into the new structure. Finally, the arrays atoms, atom_ptrs, bonds and
  bond_ptrs must be allocated to match the size of the ones in src. You should re-use (i.e. call)
  the molmalloc function in this function.
*/
molecule *molcopy(molecule *src){
  molecule *newCopy = molmalloc(src->atom_max, src->bond_max);
  for(int i = 0; i < newCopy->atom_max; i++){
    molappend_atom(newCopy, &src->atoms[i]);
  }
  for(int i = 0; i < newCopy->bond_max; i++){
    molappend_bond(newCopy, &src->bonds[i]);
  }
  //Returning the new copy
  return newCopy;
}

/*This function should free the memory associated with the molecule pointed to by ptr. This
includes the arrays atoms, atom_ptrs, bonds, bond_ptrs.*/
void molfree(molecule *ptr){
  free(ptr->atoms);
  free(ptr->atom_ptrs);
  free(ptr->bonds);
  free(ptr->bond_ptrs);
  free(ptr);
}

/*
  This function should copy the data pointed to by atom to the first “empty” atom in atoms in the
  molecule pointed to by molecule, and set the first “empty” pointer in atom_ptrs to the same
  atom in the atoms array incrementing the value of atom_no. If atom_no equals atom_max, then
  atom_max must be incremented, and the capacity of the atoms, and atom_ptrs arrays
  increased accordingly. If atom_max was 0, it should be incremented to 1, otherwise it should be
  doubled. Increasing the capacity of atoms, and atom_ptrs should be done using realloc so
  that a larger amount of memory is allocated and the existing data is copied to the new location.
*/
void molappend_atom(molecule *molecule, atom *atom){
  //if loop to check for molecule atoms and atom max
  if(molecule->atom_no == molecule->atom_max){
    if(molecule->atom_max == 0){
      molecule->atom_max++;
    }else{
      //run if condition above isnt met
      molecule->atom_max *= 2;
    }
  }
  atom_realloc(molecule);
  (molecule->atoms)[molecule->atom_no] = *atom;

  for(int i = 0; i < molecule->atom_no; i++){
    molecule->atom_ptrs[i] = &(molecule->atoms[i]);
  }
  molecule->atom_ptrs[molecule->atom_no] = &molecule->atoms[molecule->atom_no];
  //Increment molecule->atom_no++;
  molecule->atom_no++;
}

/*
This function should operate like that molappend_atom function, except for bonds.
*/
void molappend_bond(molecule *molecule, bond *bond){

  if(molecule->bond_no == molecule->bond_max){
    if(molecule->bond_max == 0){
      molecule->bond_max++;
    }else{
      //do if condition above isnt met
      molecule->bond_max *= 2;
    }
  }
  bond_realloc(molecule);
  (molecule->bonds)[molecule->bond_no] = *bond;
  for(int i = 0; i < molecule->bond_no; i++){
    molecule->bond_ptrs[i] = &(molecule->bonds[i]);
  }
  molecule->bond_ptrs[molecule->bond_no] = &molecule->bonds[molecule->bond_no];
  //increment molecule->bond_no++
  molecule->bond_no++;
}

int compare_atom_z(const void* a1, const void* a2) 
{
   //type casting for both
    atom ** atm1 = (atom **)a1;
    atom ** atm2 = (atom **)a2;
    //returning atm1 - atm2 calculation
    return (*atm1)->z - (*atm2)->z;
}

int bond_comp(const void* b1, const void* b2)
{
    //set double pointers for bnd1 and bnd2 
    bond ** bnd1 = (bond **)b1;
    bond ** bnd2 = (bond **)b2;
    //declare and set b1_avg and b2_avg and return b1_avg - b2_avg
    double b1_avg = (*bnd1)->z;
    double b2_avg = (*bnd2)->z;
    return b1_avg - b2_avg;
} 

void molsort(molecule *molecule){
  //using qsort 
  qsort(molecule->atom_ptrs, molecule->atom_no, sizeof(atom *), compare_atom_z);
  qsort(molecule->bond_ptrs, molecule->bond_no, sizeof(bond *), bond_comp);
}

void xrotation( xform_matrix xform_matrix, unsigned short deg ){
  //converting deg to rad
  long double rad = (M_PI/180) * deg;
  //calculation for matrix 
  xform_matrix[0][0] = 1;
  xform_matrix[0][1] = 0;
  xform_matrix[0][2] = 0;

  //calculation for matrix
  xform_matrix[1][0] = 0;
  xform_matrix[1][1] = cos(rad);
  xform_matrix[1][2] = -sin(rad);

  //calculation for matrix
  xform_matrix[2][0] = 0;
  xform_matrix[2][1] = sin(rad);
  xform_matrix[2][2] = cos(rad);
}

void yrotation( xform_matrix xform_matrix, unsigned short deg ){
  //converting deg to rad
  long double rad = (M_PI/180) * deg;
  //calculaiton for matrix
  xform_matrix[0][0] = cos(rad);
  xform_matrix[0][1] = 0;
  xform_matrix[0][2] = sin(rad);

  //calculation for matrix
  xform_matrix[1][0] = 0;
  xform_matrix[1][1] = 1;
  xform_matrix[1][2] = 0;

  //calculation for matrix
  xform_matrix[2][0] = -sin(rad);
  xform_matrix[2][1] = 0;
  xform_matrix[2][2] = cos(rad);
}

void zrotation( xform_matrix xform_matrix, unsigned short deg ){
  //converting deg to rad
  long double rad = (M_PI/180) * deg;
  //calculating for matrix
  xform_matrix[0][0] = cos(rad);
  xform_matrix[0][1] = -sin(rad);
  xform_matrix[0][2] = 0;

  //calculating for matrix
  xform_matrix[1][0] = sin(rad);
  xform_matrix[1][1] = cos(rad);
  xform_matrix[1][2] = 0;

  //calculating for matrix
  xform_matrix[2][0] = 0;
  xform_matrix[2][1] = 0;
  xform_matrix[2][2] = 1;
}

void mol_xform( molecule *molecule, xform_matrix matrix ){
  //declare and initialize x,y,z
  double x, y, z;
  for(int atom_i = 0; atom_i < molecule->atom_no; atom_i++){
    for(int i = 0; i < 3; i++){
      //decalre sum to 0
      double sum = 0;
      for(int j = 0; j < 3; j++){
        if(j == 0){
          sum += molecule->atoms[atom_i].x * matrix[i][j];
        }else if(j == 1){
          sum += molecule->atoms[atom_i].y * matrix[i][j];
        }else{
          sum += molecule->atoms[atom_i].z * matrix[i][j];
        }
      }
      if(i == 0){
          x = sum;
        }else if(i == 1){
          y = sum;
        }else{
          z = sum;
        }
    }
    molecule->atoms[atom_i].x = x;
    molecule->atoms[atom_i].y = y;
    molecule->atoms[atom_i].z = z;
  }
  for(int i = 0; i < molecule->bond_no; i++){
    compute_coords(&(molecule->bonds[i]));
  }
}
