/*

Copyright (C) 2003-2006 Christian Kreibich <christian@whoop.org>.
Copyright (C) 2010 Yuri Pirola <yuri.pirola (*) gmail.com>.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies of the Software and its documentation and acknowledgment shall be
given in the documentation and software packages that this Software was
used.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/
#ifndef __lst_structs_h
#define __lst_structs_h

#include <sys/queue.h>
#include "lst_string.h"
#include <stdbool.h>

#define LST_STRING_HASH_SIZE         199

typedef struct lst_stree             LST_STree;
typedef struct lst_node              LST_Node;
typedef struct lst_edge              LST_Edge;
typedef struct lst_string_hash_item  LST_StringHashItem;
typedef struct lst_string_hash       LST_StringHash;

struct lst_edge
{
  LIST_ENTRY(lst_edge)        siblings;

  LST_Node                   *src_node;
  LST_Node                   *dst_node;

  LST_StringIndex             range;
};


struct lst_node
{
  /* Each node maintains a list for its children. */
  LIST_HEAD(elist, lst_edge)  kids;

  LST_Edge                   *up_edge;

  LST_Node                   *suffix_link_node;

  int                         index;

// Campi aggiuntivi da associare al nodo
  char single_char; // Set to \0 if occurrences of subtree are preceeded
						  // by more than one symbol
  size_t* slices;
  size_t string_depth;
};


struct lst_phase_num
{
  LIST_ENTRY(lst_phase_num)   items;
  
  u_int                       phase;
};

struct lst_string_hash_item
{
  LIST_ENTRY(lst_string_hash_item) items;

  LST_String   *string;
  int           index;
};

LIST_HEAD(lst_string_hash, lst_string_hash_item);

struct lst_stree
{
  /* Number of strings currently in tree.
   */
  u_int                             num_strings;

  /* Current phase of Ukkonen's algorithm.
   * In order to implement the "Once a leaf, always
   * a leaf" Trick as explained by Gusfield, we make
   * this a pointer to an integer.
   */
  u_int                            *phase;

  /* To avoid the O(m) cost of setting that value in
   * stone once a string insertion is over, we make phase
   * point into the following list. A new element is created
   * for every string insertion.
   */
  LIST_HEAD(phase_s, lst_phase_num) phases;

  /* Current extension of Ukkonen's algorithm. */
  u_int                             ext;

  /* Well ... guess :) */
  LST_Node                         *root_node;

  /* A simple hashtable for the strings in the tree, mapping
   * them to indices starting from 1:
   */
  LST_StringHash                   *string_hash;

  /* A counter for string index numbers */
  int                               string_index;

  /* Whether or not we allow duplicates in our tree */
  int                               allow_duplicates;

  size_t* arr_occs;
};


#endif

