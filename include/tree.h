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
#ifndef _TREE_H_
#define _TREE_H_

#include "generic.h"
#include "list.h"

typedef struct _tree* ptree;
typedef struct _tree_node* ptree_node;

typedef void (*tree_visit_function)(ptree_node, item);

ptree tree_create(item root_data);

void tree_destroy(ptree tree,
						delete_function data_destroy);

ptree_node tree_get_root(ptree tree);

plist tree_node_get_children(ptree_node node);

item tree_node_get_data(ptree_node node);

void tree_node_add_child(ptree_node parent_node, item data);

ptree_node tree_node_get_parent(ptree_node node);

int tree_node_get_depth(ptree_node node);

void tree_depth_first_visit(ptree tree,
									 item aux_data,
									 tree_visit_function on_enter,
									 tree_visit_function on_exit);

void tree_print(ptree tree, FILE* f, print_function print);

#endif
