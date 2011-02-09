/**
 *
 *
 *                              PIntron
 *
 * A novel pipeline for computational gene-structure prediction based on
 * spliced alignment of expressed sequences (ESTs and mRNAs).
 *
 * Copyright (C) 2010  Yuri Pirola
 *
 * Distributed under the terms of the GNU Affero General Public License (AGPL)
 *
 *
 * This file is part of PIntron.
 *
 * PIntron is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIntron is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with PIntron.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
#include "tree.h"
#include "util.h"
#include "log.h"

struct _tree {
  ptree_node root;
};

struct _tree_node {
  void* data;
  ptree_node parent;
  plist children;
  int depth;
};

struct _tree_node_destroy {
  delete_function data_destroy;
};

static ptree_node
tree_node_create(item data);

static void
tree_node_destroy(ptree_node node,
						struct _tree_node_destroy* data_destroy);

ptree tree_create(item root_data) {
  my_assert(root_data != NULL);
  ptree tree = PALLOC(struct _tree);
  ptree_node root = tree_node_create(root_data);
  tree->root= root;
  tree->root->depth = 0;
  tree->root->parent = NULL;
  return tree;
}

void tree_destroy(ptree tree, delete_function data_destroy) {
  my_assert(tree != NULL);
  my_assert(data_destroy != NULL);
  struct _tree_node_destroy* tnd=
	 PALLOC(struct _tree_node_destroy);
  tnd->data_destroy= data_destroy;
  tree_depth_first_visit(tree,
								 (item)tnd,
								 NULL,
								 (tree_visit_function)tree_node_destroy);
  pfree(tree);
  pfree(tnd);
}

static ptree_node
tree_node_create(item data) {
  my_assert(data != NULL);
  ptree_node node = PALLOC(struct _tree_node);
  node->children = list_create();
  node->data= data;
  node->depth= 0;
  node->parent= NULL;
  return node;
}

void
tree_node_destroy(ptree_node node,
						struct _tree_node_destroy* tnd) {
  my_assert(node != NULL);
  my_assert(tnd != NULL);
  tnd->data_destroy(node->data);
  list_destroy(node->children, noop_free);
  pfree(node);
}

ptree_node tree_get_root(ptree tree) {
  my_assert(tree != NULL);
  return tree->root;
}

plist tree_node_get_children(ptree_node node) {
  my_assert(node != NULL);
  return node->children;
}

item tree_node_get_data(ptree_node node) {
  my_assert(node != NULL);
  return node->data;
}

void tree_node_add_child(ptree_node parent_node, item data) {
  my_assert(parent_node != NULL);
  my_assert(data != NULL);
  ptree_node new_node= tree_node_create(data);
  new_node->parent= parent_node;
  new_node->depth= (parent_node->depth) + 1;
  list_add_to_tail(parent_node->children, new_node);
}

int tree_node_get_depth(ptree_node node) {
  my_assert(node != NULL);
  return node->depth;
}

static void
tree_node_depth_first_visit(ptree_node node,
									 item aux_data,
									 tree_visit_function on_enter,
									 tree_visit_function on_exit) {

  my_assert(node != NULL);

  if (on_enter!=NULL)
	 on_enter(node, aux_data);

  plistit it = list_first(tree_node_get_children(node));
  while (listit_has_next(it)){
	 tree_node_depth_first_visit((ptree_node)listit_next(it),
										  aux_data,
										  on_enter,
										  on_exit);
  }
  listit_destroy(it);
  if (on_exit!=NULL)
	 on_exit(node, aux_data);
}

void
tree_depth_first_visit(ptree tree,
							  item aux_data,
							  tree_visit_function on_enter,
							  tree_visit_function on_exit) {

  my_assert(tree != NULL);

  tree_node_depth_first_visit(tree_get_root(tree),
										aux_data,
										on_enter,
										on_exit);
}

struct _print_info {
  print_function print;
  FILE* f;
};

static void
tree_node_print_onenter(ptree_node node, struct _print_info* pi) {

  my_assert(node != NULL);
  my_assert(pi != NULL);

  print_repetitions(pi->f, ' ', tree_node_get_depth(node));
  fprintf(pi->f, "(");
  pi->print(pi->f, tree_node_get_data(node));
  fprintf(pi->f, "\n");
}

static void
tree_node_print_onexit(ptree_node node, struct _print_info* pi) {

  my_assert(node != NULL);
  my_assert(pi != NULL);

  print_repetitions(pi->f, ' ', tree_node_get_depth(node));
  fprintf(pi->f, ")\n");
}

void tree_print(ptree tree, FILE* f, print_function print) {

  my_assert(tree != NULL);
  my_assert(f != NULL);
  my_assert(print != NULL);

  struct _print_info* pi= PALLOC(struct _print_info);
  pi->f= f;
  pi->print= print;
  tree_depth_first_visit(tree, pi,
								 (tree_visit_function)tree_node_print_onenter,
								 (tree_visit_function)tree_node_print_onexit);
  pfree(pi);
}
