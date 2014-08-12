/* libavl - manipulates AVL trees.
   Copyright (C) 1998, 1999 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   02111-1307, USA.

   The author may be contacted at <pfaffben@pilot.msu.edu> on the
   Internet, or as Ben Pfaff, 12167 Airport Rd, DeWitt MI 48820, USA
   through more mundane means. */

/* This is file avltr.h in libavl. */

/* Martin Schlather: 21 October, 2011,
   This file has been file avltr.h in libavl. 
   The last 12 lines have been put into / *  * /
   as they are incompatible with R.
   The above stated GNU General Public License to its full extend
   is still valid.
*/

#include"RF.h"

#if !avltr_h
#define avltr_h 1

/* The default maximum height of 32 allows for AVL trees having
   between 5,704,880 and 4,294,967,295 nodes, depending on order of
   insertion.  You may change this compile-time constant as you
   wish. */
#ifndef AVL_MAX_HEIGHT
#define AVL_MAX_HEIGHT	32
#endif


/* Structure for a node in a right-threaded AVL tree. */
typedef struct avltr_node
  {
    cell_type  *data;			/* Pointer to data. */
    struct avltr_node *link[2];	/* Subtrees or threads. */
    signed char bal;		/* Balance factor. */
    char cache;			/* Used during insertion. */
    char pad;		        /* Reserved for fully threaded trees. */
    signed char rtag;		/* Right thread tag. */
  }
avltr_node;

/* Used for traversing a right-threaded AVL tree. */
typedef struct avltr_traverser
  {
    int init;				/* Initialized? */
    avltr_node *p;		/* Last node returned. */
  }
avltr_traverser;

/* Initializer for avltr_traverser. */
#define AVLTR_TRAVERSER_INIT {0}

/* Function types. */
#if !AVL_FUNC_TYPES
#define AVL_FUNC_TYPES 1
typedef int (*avl_comparison_func) (const void *a, const void *b, void *param);
typedef void (*avl_node_func) (void *data, void *param);
typedef void *(*avl_copy_func) (void *data, void *param);
#endif

/* Structure which holds information about a threaded AVL tree. */
typedef struct avltr_tree
  {
    avltr_node root;		/* Tree root node. */
    avl_comparison_func cmp;	/* Used to compare keys. */
    int count;			/* Number of nodes in the tree. */
    int *param;		/* Arbitary user data. */
  }
avltr_tree;

/* General functions. */
 avltr_tree *avltr_create (avl_comparison_func, int *param);
 void avltr_destroy (avltr_tree *, avl_node_func);
 void avltr_free (avltr_tree *);
 int avltr_count (avltr_tree *);
 avltr_tree *avltr_copy (avltr_tree *, avl_copy_func);
struct avl_tree;
 avltr_tree *avltr_thread (struct avl_tree *);
 struct avl_tree *avltr_unthread (avltr_tree *);

/* Walk the tree. */
 void avltr_walk (avltr_tree *, avl_node_func, int *param);
 cell_type *avltr_traverse (avltr_tree *, avltr_traverser *);
#define avlt_init_traverser(TRAVERSER) ((TRAVERSER)->init = 0)
 cell_type **avltr_next (avltr_tree *tree, cell_type **item);

/* Search for a given item. */
 cell_type **avltr_probe (avltr_tree *, cell_type *);
 cell_type *avltr_delete (avltr_tree *, cell_type *);
 cell_type **avltr_find (avltr_tree *, cell_type *);
 cell_type **avltr_find_close (avltr_tree *, cell_type *);

#if __GCC__ >= 2
 inline cell_type *avltr_insert (avltr_tree *tree, cell_type *item)
{
  cell_type **p = avltr_probe (tree, item);
  return (*p == item) ? NULL : *p;
}

 inline cell_type *avltr_replace (avltr_tree *tree, cell_type *item)
{
  cell_type **p = avltr_probe (tree, item);
  if (*p == item)
    return NULL;
  else
    {
      cell_type *r = *p;
      *p = item;
      return r;
    }
}
#else /* not gcc */
 cell_type *avltr_insert (avltr_tree *tree,  cell_type *item);
 cell_type *avltr_replace (avltr_tree *tree, cell_type *item);
#endif /* not gcc */

/* Easy assertions on insertion & deletion. */
/* *** Martin Schlather: uncomment the following lines
   to get the original version back:

#ifndef NDEBUG
#define avltr_force_insert(A, B)		\
	do					\
	  {					\
            void *r = avltr_insert (A, B);	\
	    assert (r == NULL);			\
	  }					\
	while (0)
 void *avltr_force_delete (avltr_tree *, void *);
#else
#define avltr_force_insert(A, B)		\
	avltr_insert (A, B)
#define avltr_force_delete(A, B)		\
	avltr_delete (A, B)
#endif

*/

#endif /* avltr_h */


