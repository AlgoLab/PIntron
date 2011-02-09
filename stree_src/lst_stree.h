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
#ifndef __lst_stree_h
#define __lst_stree_h

#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/queue.h>

#include "lst_structs.h"


/**
 * lst_stree_new - creates a suffix tree for a set of strings.
 * @strings: set of strings to build tree with.
 *
 * This is an implementation of Ukkonen's O(n) algorithm for creating
 * a suffix tree. Upon return, the tree contains information on all
 * the strings contained in the given string set. If you don't want
 * to insert strings right away, just pass %NULL.
 *
 * Returns: new suffix tree.
 */
LST_STree   *lst_stree_new(LST_StringSet *strings);


/**
 * lst_stree_free - suffix tree destructor.
 * @tree: tree to clean up.
 *
 * The function releases all the memory claimed by the suffix tree.
 * It does not touch any of the strings contained in the tree when
 * called, it only cleans up the tree itself. Use when the tree
 * was created with lst_stree_new().
 */
void         lst_stree_free(LST_STree *tree);


/**
 * lst_stree_init - suffix tree initialization for existing tree structure.
 * @tree: tree structure to initialize.
 *
 * This function initializes a tree structure that already exists.
 * It is hence faster when you need a suffix tree in a tight loop
 * as no data need be allocated and later on freed. It does not check
 * if any data is existing in the structure when called; make sure you
 * call lst_stree_clear() when you want to use the structure repeatedly.
 *
 * Returns: value > 0 when initialization was successful, 0 otherwise.
 */
int          lst_stree_init(LST_STree *tree);


/**
 * lst_stree_clear - cleans up internal tree structure.
 * @tree: tree to clear.
 *
 * This is the counterpart to lst_stree_init(). It cleans up the tree
 * but does not free the tree structure itself.
 */
void         lst_stree_clear(LST_STree *tree);


/**
 * lst_stree_add_string - adds a string from to tree.
 * @tree: tree to add string to.
 * @string: string to add.
 *
 * The function adds @string to the tree, unless he string is
 * a duplicate of an existing string and duplicates are not
 * allowed (see lst_stree_allow_duplicates()).
 */
void         lst_stree_add_string(LST_STree *tree, LST_String *string);


/**
 * lst_stree_remove_string - removes a string from the tree.
 * @tree: tree to remove string from.
 * @string: string to remove.
 *
 * The function checks whether @tree in fact contains @string and
 * if that's the case, removes it from the tree.
 */
void         lst_stree_remove_string(LST_STree *tree, LST_String *string);


/**
 * lst_stree_get_string_index - returns a nonnegative index for a string.
 * @tree: tree to query.
 * @string: string to look up.
 *
 * Within a suffix tree, every string contained in it is associated with
 * an integer index value. This function returns that value.
 *
 * Returns: index of @string in @tree.
 */
int          lst_stree_get_string_index(LST_STree *tree, LST_String *string);


/**
 * lst_stree_allow_duplicates - whether the tree may contain duplicates.
 * @tree: tree to modify.
 * @duplicates_flag: whether to allow duplicates (> 0) or not (0).
 *
 * Depending on the application of the suffix tree, it may be okay to
 * have duplicates of strings in the tree or not. By default, duplicates
 * are allowed. However, if you want to prevent insertion of a string
 * that is already contained in the tree, pass 0.
 */
void         lst_stree_allow_duplicates(LST_STree *tree, int duplicates_flag);


/**
 * lst_node_get_parent - returns parent of a node.
 * @node: node to find parent for.
 *
 * Returns: the parent node of a node, or %NULL if no
 * such node exists.
 */
LST_Node    *lst_node_get_parent(LST_Node *node);


/**
 * lst_node_is_leaf - checks whether a node is a leaf.
 * @node: node to check.
 *
 * Returns: value > 0 if @node is a leaf, 0 otherwise.
 */
int lst_node_is_leaf(LST_Node *node);


/**
 * lst_node_is_root - checks whether a node is the tree root.
 * @node: node to check.
 *
 * Returns: value > 0 if @node is the root, 0 otherwise.
 */
int lst_node_is_root(LST_Node *node);


/**
 * lst_node_get_string_length - returns length of string leading to node.
 * @node: node to query.
 *
 * Returns: the number of string items found on the edges iterated
 * when going from the root down to @node.
 */
int lst_node_get_string_length(LST_Node *node);


/**
 * lst_edge_get_length - returns the length of a substring on an edge.
 * @edge: edge to query.
 *
 * Returns: the length of the substring associated with that edge.
 */
size_t lst_edge_get_length(const LST_Edge * const edge);

#endif
