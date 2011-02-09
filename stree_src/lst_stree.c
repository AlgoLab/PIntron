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
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/queue.h>

#include "lst_stree.h"
#include "lst_string.h"


extern int string_id_counter;


/* The number of strings for which we initially
 * make room in the slots array
 */
#define LST_STREE_STRINGSLOTS   32


/* A path in an implicit suffix tree can end at either a node, or
 * at some point in the label of an edge. We remember in each
 * extension where we ended using an LST_PathEnd structure.
 */
typedef struct lst_path_end
{
  LST_Node    *node;

  LST_Edge    *edge;
  u_int        offset;

} LST_PathEnd;


static LST_Edge *
edge_new(LST_Node *src_node,
	 LST_Node *dst_node,
	 LST_String *string,
	 u_int start_index, u_int end_index)
{
  LST_Edge *edge;
  
  edge = malloc(sizeof(LST_Edge));

  if (!edge)
    {
      return NULL;
    }

  edge->src_node = src_node;
  edge->dst_node = dst_node;

  edge->range.string = string;
  edge->range.start_index = start_index;
  edge->range.end_index_local = end_index;
  edge->range.end_index = &edge->range.end_index_local;
  
  dst_node->up_edge = edge;

  LIST_INSERT_HEAD(&src_node->kids, edge, siblings);

  return edge;
}


static LST_Edge *
edge_leaf_new(LST_STree *tree,
	      LST_Node *src_node,
	      LST_Node *dst_node,
	      LST_String *string,
	      u_int start_index)
{
  LST_Edge *edge = edge_new(src_node, dst_node, string, start_index, 0);

  if (!edge)
    return NULL;

  /* For a leaf edge, we make the end index point to
   * the current phase pointer in the tree structure.
   */
  edge->range.end_index = tree->phase;
  
  return edge;
}


static void
edge_free(LST_Edge *edge)
{
  if (edge)
    free(edge);
}


static LST_Node *
node_new(int index)
{
  LST_Node *node;

  node = malloc(sizeof(LST_Node));

  if (!node)
    {
      return NULL;
    }

  node->up_edge= NULL;
  node->suffix_link_node= NULL;
  LIST_INIT(&node->kids);
  node->index = index;

  return node;
}


static void
node_free(LST_Node *node)
{
  LST_Edge *edge;

  if (!node)
    return;

  /* Recursively clean up kids first */
  while (node->kids.lh_first)
    {
      node_free(node->kids.lh_first->dst_node);      
      edge = node->kids.lh_first;
      LIST_REMOVE(node->kids.lh_first, siblings);
      edge_free(edge);
    }

  if (node->slices!=NULL)
	 free(node->slices);
  free(node);
}

static LST_Edge *
node_find_edge_with_startitem(LST_Node *node, LST_String *string, u_int index)
{
  LST_Edge *edge = NULL;

  if (!node || !string || index >= string->num_items)
    {
      return NULL;
    }

  for (edge = node->kids.lh_first; edge; edge = edge->siblings.le_next)
    {
      /* Skip this edge if the first characters don't match
       * what we're looking for.
       */
      if (lst_string_eq(edge->range.string, edge->range.start_index,
			string, index))
	return edge;
    }

  return NULL;
}


static void
stree_extend_label_at_leaf(LST_Node *leaf)
{
  if (!lst_node_is_leaf(leaf))
    return;
  
  *leaf->up_edge->range.end_index = *leaf->up_edge->range.end_index + 1;
}


static void
stree_path_end_node(LST_PathEnd *end, LST_Node *node)
{
  memset(end, 0, sizeof(LST_PathEnd));
  end->node = node;
}


static void
stree_path_end_edge(LST_PathEnd *end, LST_Edge *edge, u_int offset)
{
  memset(end, 0, sizeof(LST_PathEnd));
  end->edge   = edge;
  end->offset = offset;
}


static void
stree_path_end_advance(LST_PathEnd *end, LST_Edge *edge)
{
  if (!end)
    return;

  if (end->node)
    {
      if (lst_edge_get_length(edge) == 1)
	stree_path_end_node(end, edge->dst_node);
      else
	stree_path_end_edge(end, edge, 1);
    }
  else
    {
      end->offset++;
      
      if (end->offset == (u_int) lst_edge_get_length(end->edge))
	stree_path_end_node(end, end->edge->dst_node);
    }
}


/**
 * stree_follow_string_slow - follows an arbitrary string through the tree.
 * @tree: tree to query.
 * @node: node to start at.
 * @string: the string to follow.
 * @end: result argument, see below.
 *
 * The function follows the @skipstring in @tree starting from @node,
 * reporting where in the tree the string ends through the @end
 * result pointer.
 *
 * Returns: number of string items successfully matched.
 */
static u_int
stree_follow_string_slow(LST_STree *tree, LST_Node *node,
			 LST_String *string, LST_PathEnd *end)
{
  LST_Edge   *edge = NULL;
  u_int       items_todo = 0, items_done = 0, common;

  if (!tree || !node || !string || !end)
    {
      memset(end, 0, sizeof(LST_PathEnd));
      return 0;
    }
  
  items_todo = string->num_items;
      
  /* Find edge where our string starts, making use of the fact
   * that no two out-edges from a node can start with the same
   * character.
   */
  while (items_todo > 0)
    {
      edge = node_find_edge_with_startitem(node, string, items_done);

      if (!edge)
	{
	  stree_path_end_node(end, node);
	  return items_done;
	}

      common =
	lst_string_items_common(edge->range.string, edge->range.start_index,
				string, items_done,
				items_todo);

      if (common < (u_int) lst_edge_get_length(edge))
	{
	  stree_path_end_edge(end, edge, common);
	  return items_done + common;
	}
      
		node = edge->dst_node;
		items_done += lst_edge_get_length(edge);
		items_todo -= lst_edge_get_length(edge);
    }
  
  stree_path_end_node(end, node);
  return items_done;
}


/**
 * stree_follow_string - follows an existing string in the tree, using skip/count.
 * @tree: tree to query.
 * @node: node to start at.
 * @skipstring: the string to follow.
 * @end: result argument, see below.
 *
 * The function follows the @skipstring in @tree starting from @node,
 * reporting where in the tree the string ends through the @end
 * result pointer.
 */
static void
stree_follow_string(LST_STree *tree, LST_Node *node,
						  LST_StringIndex *skipstring, LST_PathEnd *end)
{
  LST_Edge   *edge = NULL;
  LST_String *string;
  u_int       items_todo = 0, items_done = 0, common;

  if (!tree || !node || !skipstring || !end) {
	 return;
  }

  string = skipstring->string;

/* We need to figure out how many string items we need to walk down in
   * the tree so that we can then extend by the next item.
   */
  if (skipstring->start_index == LST_EMPTY_STRING) {
	 stree_path_end_node(end, node);
	 return;
  }

  items_todo = *(skipstring->end_index) - skipstring->start_index + 1;

  if (items_todo == 0) {
	 stree_path_end_node(end, node);
	 return;
  }

  /* Find edge where our string starts, making use of the fact
   * that no two out-edges from a node can start with the same
   * character.
   */
  while (items_todo > 0)
    {
      u_int edge_len;

      edge = node_find_edge_with_startitem(node, string, skipstring->start_index + items_done);
      
      if (!edge)
	{
	  stree_path_end_node(end, node);
	  return;
	}

      edge_len = (u_int) lst_edge_get_length(edge);

      /* Follow edges in tree, emplying the Skip/Count Trick as per Gusfield.
       * When the string we're looking up is longer than the edge's label,
       * we can just skip this edge:
       */
		if (items_todo >= edge_len) {
		  items_todo -= edge_len;
		  items_done += edge_len;

		  if (items_todo == 0) {
			 stree_path_end_node(end, edge->dst_node);
			 return;
		  }

		  node = edge->dst_node;
		  continue;
		}

/* When the string is shorter than the edge, we need to compare
 * the strings and figure out where we stop within the edge.
 * This will need a new edge as per extension Rule 2.
 */
		break;
	 }

  common =
	 lst_string_items_common(edge->range.string,
									 edge->range.start_index,
									 string,
									 skipstring->start_index + items_done,
									 items_todo);
  stree_path_end_edge(end, edge, common);
}


/**
 * stree_get_skipstring - finds the string range we need to cross up to previous node.
 * @tree: tree in which we search.
 * @end: definition of the tree location where we start going up.
 * @skipstring: result pointer, see below.
 *
 * The function finds the string range we need to jump over until
 * we get back up to a node that has a suffix link. If that node is
 * the root, the upstring will be empty.
 *
 * Returns: The node we arrive at. Also, the string range is returned
 * through @skipstring.
 */
static LST_Node *
stree_get_skipstring(LST_STree *tree, LST_PathEnd *end, LST_StringIndex *skipstring)
{
  LST_Node *parent_node;

  if (end->node)
    {
      if (end->node->suffix_link_node)
	{
	  /* The node we ended at already has a suffix link,
	   * so we need to do nothing. Mark the skipstring
	   * as empty:
	   */
	  skipstring->start_index = LST_EMPTY_STRING;

	  return end->node->suffix_link_node;
	}

      /* If the node doesn't have a suffix link directly, we must
       * hop up over the complete edge's string. If we end up at
       * the root, the caller must descend as in the naive algorithm,
       * otherwise the skipstring is just whatever the string attached
       * to the edge is:
       */
      if (end->node == tree->root_node)
	{
	  skipstring->start_index = LST_EMPTY_STRING;
	  return tree->root_node;
	}

      parent_node = lst_node_get_parent(end->node);
      
      if (parent_node == tree->root_node)
	{
	  return tree->root_node;
	}

      lst_string_index_copy(&(end->node->up_edge->range), skipstring);
      
      /* Follow the node up the edge, and cross over across
       * the suffix link. Then return the node we arrive at.
       */
      return parent_node->suffix_link_node;
    }

  /* Okay -- the funkier case: we start in the middle of an edge. */
  parent_node = end->edge->src_node;
  
  if (parent_node == tree->root_node)
    {
      return tree->root_node;
    }
  
  skipstring->string = end->edge->range.string;
  skipstring->start_index = end->edge->range.start_index;
  *(skipstring->end_index) = skipstring->start_index + end->offset - 1;
  

  return parent_node->suffix_link_node;
}


/**
 * stree_extend_at_node - extends a tree by one string item, from a node.
 * @tree: tree to extend.
 * @string: string containing the new item.
 * @end: definition of where to extend from.
 * @extend_index: index in @string to insert into tree.
 * @stop_extensions: result pointer, see below.
 *
 * The function extends the tree by a single string item, down from
 * a node within the tree. The point of insertion is given through @end.
 *
 * Returns: 1 in @stop_extensions if Rule 3 (see Gusfield) applied and we
 * can hence stop extensions for the current phase.
 */
static void
stree_extend_at_node(LST_STree *tree, LST_String *string,
		     LST_PathEnd *end, u_int extend_index, u_int *stop_extensions)
{
  LST_Edge *edge;
  LST_Node *new_node;

  if (lst_node_is_leaf(end->node))
    {
      stree_extend_label_at_leaf(end->node);
      stree_path_end_edge(end, end->node->up_edge,
			  lst_edge_get_length(end->node->up_edge) - 2);
    }
  else
    {
      edge = node_find_edge_with_startitem(end->node, string, extend_index);
      
      if (!edge)
	{
	  /* Extension Rule 2: */
	  new_node = node_new(tree->ext);
	  edge = edge_leaf_new(tree, end->node, new_node,
			       string, extend_index);
	}
      else
	{
	  /* Otherwise it's Extension Rule 3, so we do nothing,
	   * but only mark that we've applied that rule so we
	   * can speed things up. */
	  stree_path_end_advance(end, edge);
	  *stop_extensions = 1;
	}
    }
}


/**
 * stree_extend_at_edge - extends a tree by one string item, from within an edge.
 * @tree: tree to extend.
 * @string: string containing the new item.
 * @end: definition of where to extend from.
 * @extend_index: index in @string to insert into tree.
 * @stop_extensions: result pointer, see below.
 *
 * The function extends the tree by a single string item, down from within
 * an edge in the tree. The point of insertion is given through @end.
 *
 * Returns: 1 in @stop_extensions if Rule 3 (see Gusfield) applied and we
 * can hence stop extensions for the current phase.
 */
static LST_Node *
stree_extend_at_edge(LST_STree *tree, LST_String *string,
							LST_PathEnd *end, u_int extend_index, u_int *stop_extensions)
{
  LST_Node *old_node, *new_inner_node = NULL, *new_node;
  LST_Edge *new_edge;

  if (lst_string_eq(string, extend_index,
						  end->edge->range.string, end->edge->range.start_index + end->offset))
	 {
		stree_path_end_advance(end, NULL);
		*stop_extensions = 1;
	 }
  else
	 {
		new_inner_node = node_new(-1);
		old_node = end->edge->dst_node;
/* Carefully carefully carefully -- when we split a leaf edge,
 * we need to figure out what kind of end index to use (the edge-local
 * or global one). It's not enough to check whether it's a leaf or
 * not -- it could be a leaf created for a previous string. So only
 * make the end index point back into the tree if it was pointing there
 * originally anyway. However, tree->phase changes over time, so we
 * must not use that.
 */
		if (lst_node_is_leaf(old_node))
		  {
			 u_int *end_index = end->edge->range.end_index;

			 new_edge =
				edge_leaf_new(tree, new_inner_node, old_node,
								  end->edge->range.string,
								  end->edge->range.start_index + end->offset);

			 new_edge->range.end_index = end_index;
		  }
		else
		  {
			 new_edge =
				edge_new(new_inner_node, old_node,
							end->edge->range.string,
							end->edge->range.start_index + end->offset,
							*(end->edge->range.end_index));
		  }

		end->edge->range.end_index = &(end->edge->range.end_index_local);
		*(end->edge->range.end_index) = end->edge->range.start_index + end->offset - 1;
		end->edge->dst_node = new_inner_node;
		new_inner_node->up_edge = end->edge;

/* Now add another edge to the new node inserted, and
 * label it with the remainder of the string.
 */
		new_node = node_new(tree->ext);
		new_edge = edge_leaf_new(tree, new_inner_node, new_node,
										 string, extend_index);
	 }

  return new_inner_node;
}


/**
 * stree_next_phase - Trick 3 management.
 * @tree: tree to advance phase pointer in.
 *
 * This function takes care of allowing Trick 3 to work correctly
 * when multiple strings are inserted into the tree, and helps avoiding
 * the O(m) overhead to create a true suffix tree from the intermediate
 * one Gusfield mentions in section 6.1.6. We don't exchange any labels
 * when we consider edge label ends "final", we rather start using
 * different pointers to end offsets from now on, not touching the old
 * ones any more in the future.
 */
static void
stree_next_phase(LST_STree *tree)
{
  struct lst_phase_num *phase;

  if (!tree)
	 return;

  phase = malloc(sizeof(struct lst_phase_num));
  LIST_INSERT_HEAD(&tree->phases, phase, items);
  tree->phase = &phase->phase;

}


/**
 * stree_add_string - inserts a new string into the tree.
 * @tree: tree to insert into.
 * @string: new string to insert.
 *
 * The function inserts a string into the tree. It recognizes whether
 * it's the first or an additional string and applies Extension Rules
 * 1-3 accordingly.
 */
static void
stree_add_string_impl(LST_STree *tree, LST_String *string)
{
  LST_PathEnd end;
  LST_Node *node, *start_leaf = NULL, *inner_node;
  LST_Node *prev_new_node = NULL, *new_inner_node;
  LST_Edge *edge;
  LST_StringIndex skipstring;
  u_int i, j, stop_extensions = 0, last_extension = 0, use_end = 0, find_skipstring = 1;

  if (!tree || !string)
	 return;

  stree_next_phase(tree);

/* First and subsequent string insertions are handled slightly
 * differently -- first one is built from scratch, the others
 * are built after matching them onto the current tree as much
 * as possible, and then proceeding as usual.
 */
  if (tree->num_strings == 0)
	 {

		node = node_new(0);
		edge = edge_leaf_new(tree, tree->root_node, node, string, 0);

/* Construct first tree -- simply enter a single edge with
 * a single-item label.
 */

		start_leaf = node;
		stree_path_end_edge(&end, start_leaf->up_edge, 0);
		last_extension = 0;

/* Skip to phase 1: */
		i = 1;
	 }
  else
	 {
/* Follow new string through tree, hence finding the number
 * of phases we can skip:
 */
      
		i = stree_follow_string_slow(tree, tree->root_node, string, &end);
		last_extension = -1;
		find_skipstring = 0;

	 }

/* Following Ukkonen's algorithm, we insert the new string
 * of length |m| in m "phases", constructing tree i+1 from
 * tree i. Each phase i+1 consists of i+1 "extensions", one
 * for each suffix in the string s[1,i+1].
 *
 * By including the last string item (our end-of-string marker),
 * we automatically avoid the issues of having suffixes that
 * are prefixes of other suffixes.
 */
  for (/* i inited above */; i < string->num_items; i++)
	 {

		*(tree->phase) = i;

		use_end = 0;
		if (stop_extensions)
		  use_end = 1;
		stop_extensions = 0;
      
      
/* Now do the remaining extensions. We don't start at index 0
 * with extensions but rather implement speedup Trick 3 as per Gusfield.
 */
		for (j = last_extension + 1; j < i + 1 && !stop_extensions; j++) {
		  tree->ext = j;
/* Get the node from which we start to walk down,
 * either found via suffix links or because it's the root:
 */
		  lst_string_index_init(&skipstring);
		  skipstring.string = string;
		  u_int extra_index = i;

		  if (use_end) {
			 use_end = 0;
		  } else if (find_skipstring) {
			 node = stree_get_skipstring(tree, &end, &skipstring);

			 if (node == tree->root_node) {
				if (skipstring.start_index != LST_EMPTY_STRING) {
/* It's the root node -- just follow the path down
 * in the tree as in the naive algorithm.
 */
				  skipstring.string = string;
				  skipstring.start_index  = j;
				  *(skipstring.end_index) = i-1;

				  if (i-1 < j)
					 skipstring.start_index = LST_EMPTY_STRING;
				}
			 } else {
/* It's not the root node -- exploit suffix
 * links as much as possible. stree_get_skipstring()
 * has already filled in stuff for us.
 */
			 }
			 stree_follow_string(tree, node, &skipstring, &end);
		  }
		  find_skipstring = 1;
		  new_inner_node = NULL;
		  inner_node = NULL;
/* Now extend the found path in the tree by the new character
 * in this phase:
 */
		  if (end.node) {
/* We followed the path up to a node. If that node is a leaf,
 * we're done per Extension Rule 1. Otherwise we need to find
 * out if our new character is somewhere on the out-edges
 * of that node, and if not, hook in a new edge.
 */
			 inner_node = end.node;
			 stree_extend_at_node(tree, string, &end,
										 extra_index,
										 &stop_extensions);
		  } else {
/* We followed the path down to somewhere within an edge
 * label. Now we need to check if the item following the
 * common part in the label is what we want to add; that's
 * a case of Extension Rule 3. Otherwise we need to split
 * the edge and hook in a new edge with the new character.
 */
			 new_inner_node=
				stree_extend_at_edge(tree, string,
											&end, extra_index,
											&stop_extensions);
		  }
	  
/* Now take care of suffix links: if we've created a new inner
 * node in the last extension, create a suffix link to either
 * the last inner node we've encountered above, or a new inner
 * node, if we've created one.
 */
		  if (prev_new_node)
			 prev_new_node->suffix_link_node = (new_inner_node ? new_inner_node : inner_node);
	  
		  prev_new_node = new_inner_node;
		}

/* Remember how far we got for the next phase. Also,
 * repeat the last extension if we aborted due to Rule 3.
 */
		last_extension = j - 1 - stop_extensions;
	 }

  tree->num_strings++;
}


void
lst_stree_add_string(LST_STree *tree, LST_String *string)
{
  LST_StringHashItem *item;
  LST_StringHash *hashlist;
  LST_PathEnd end;
  u_int items_done;

  if (!tree || !string)
	 return;

  if (! tree->allow_duplicates)
	 {
		items_done = stree_follow_string_slow(tree, tree->root_node, string, &end);
		if (items_done == string->num_items - 1)
		  return;
	 }

  item = calloc(1, sizeof(LST_StringHashItem));
  hashlist = &tree->string_hash[string->id % LST_STRING_HASH_SIZE];
  item->string = string;
  item->index  = tree->string_index++;
  LIST_INSERT_HEAD(hashlist, item, items);
  
  stree_add_string_impl(tree, string);
}


LST_STree   *
lst_stree_new(LST_StringSet *strings)
{
  LST_STree  *tree = NULL;
  LST_String *string;

  if (! (tree = malloc(sizeof(LST_STree))))
	 {
		return NULL;
	 }

  if (!lst_stree_init(tree))
	 {
		free(tree);
		return NULL;
	 }

  string_id_counter= 0;
  if (strings)
	 {
		for (string = strings->members.lh_first; string; string = string->set.le_next)
		  lst_stree_add_string(tree, string);
	 }
  
  return tree;
}


void
lst_stree_free(LST_STree *tree)
{
  if (!tree)
	 return;
  
  lst_stree_clear(tree);
  free(tree);
}


int
lst_stree_init(LST_STree *tree)
{
  int i;

  if (!tree)
	 return 0;

  memset(tree, 0, sizeof(LST_STree));

  tree->allow_duplicates = 1;
  LIST_INIT(&tree->phases);

  tree->root_node = node_new(-1);
  if (!tree->root_node)
	 goto error_return;

  tree->string_hash = calloc(LST_STRING_HASH_SIZE, sizeof(LST_StringHash));
  if (!tree->string_hash)
	 goto error_return;

  for (i = 0; i < LST_STRING_HASH_SIZE; i++)
	 LIST_INIT(&tree->string_hash[i]);

  return 1;

 error_return:

  if (tree->root_node)
	 node_free(tree->root_node);

  if (tree->string_hash)
	 free(tree->string_hash);

  return 0;
}


void
lst_stree_clear(LST_STree *tree)
{
  LST_StringHash *hash;
  LST_StringHashItem *hi;
  struct lst_phase_num *phase;
  int i;

/* Clean up the tree itself */
  node_free(tree->root_node);

/* Clean up the phases array */
  while (tree->phases.lh_first)
	 {
		phase = tree->phases.lh_first;
		LIST_REMOVE(tree->phases.lh_first, items);
		free(phase);
	 }

/* Clean up string hash */
  for (i = 0; i < LST_STRING_HASH_SIZE; i++)
	 {
		hash = &tree->string_hash[i];

		while (hash->lh_first)
		  {
			 hi = hash->lh_first;
			 LIST_REMOVE(hash->lh_first, items);
			 lst_string_free(hi->string);
			 free(hi);
		  }
	 }
  free(tree->string_hash);
}



void
lst_stree_allow_duplicates(LST_STree *tree, int duplicates_flag)
{
  if (!tree)
	 return;

  tree->allow_duplicates = duplicates_flag;
}


int
lst_stree_get_string_index(LST_STree *tree, LST_String *string)
{
  LST_StringHash *hash;
  LST_StringHashItem *hi;

  if (!tree || !string)
	 return -1;

  hash = &tree->string_hash[string->id % LST_STRING_HASH_SIZE];
  
  for (hi = hash->lh_first; hi; hi = hi->items.le_next)
	 {
		if (hi->string->id == string->id)
		  return hi->index;
	 }

  return -1;
}


LST_Node *
lst_node_get_parent(LST_Node *node)
{
  if (node->up_edge == NULL)
	 return NULL;

  return node->up_edge->src_node;
}


int
lst_node_is_leaf(LST_Node *node)
{
  if (!node)
	 return 0;

  return (node->kids.lh_first == NULL);
}


int
lst_node_is_root(LST_Node *node)
{
  if (!node)
	 return 0;

  return (node->up_edge == NULL);
}


int
lst_node_get_string_length(LST_Node *node)
{
  int depth = 0;

  if (!node)
	 return 0;

  while (! lst_node_is_root(node))
	 {
		depth += lst_edge_get_length(node->up_edge);
		node = node->up_edge->src_node;
	 }

  return depth;
}


size_t
lst_edge_get_length(const LST_Edge * const edge)
{
  if (!edge)
	 return 0;

  return *(edge->range.end_index) - edge->range.start_index + 1;
}


