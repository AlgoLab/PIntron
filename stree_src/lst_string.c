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

#include <sys/param.h>

#include "lst_string.h"
#include "aug_suffix_tree.h"

int string_id_counter;

void
lst_string_init(LST_String *string, void *data, u_int item_size, u_int num_items) {
  if (!string || !data || item_size == 0) {
	 return;
  }

  memset(string, 0, sizeof(LST_String));
  string->id = ++string_id_counter;

  string->data= data;
  string->data_external = 1;
  string->num_items= num_items + 1;
  string->item_size= item_size;
}


void *
lst_string_get_item(LST_String *string, u_int index)
{
  char *data = (char *) string->data;

  if (index >= string->num_items)
	 return NULL;

  return (void *) (data + (index * string->item_size));
}


void
lst_string_free(LST_String *string)
{
  if (!string)
	 return;

  if (string->data && !string->data_external)
	 free(string->data);

  free(string);
}


void *
lst_string_free_keep_data(LST_String *string)
{
  void *data;

  if (!string)
	 return NULL;

  data = string->data;
  free(string);

  return data;
}


u_int
lst_string_get_length(LST_String *string)
{
  if (!string)
	 return 0;

  return string->num_items - 1;
}


const char *
lst_string_print(LST_String *string)
{
  LST_StringIndex tmp_range;

  if (!string)
	 return NULL;

  lst_string_index_init(&tmp_range);

  tmp_range.string = string;
  tmp_range.start_index  = 0;
  *(tmp_range.end_index) = string->num_items - 1;

  return string_print_func_with_dollars(&tmp_range);
}


void
lst_string_item_copy(LST_String *src, u_int src_index,
							LST_String *dst, u_int dst_index)
{
  void *src_item, *dst_item;

  if (!src || !dst || src_index >= src->num_items || dst_index >= dst->num_items)
	 return;

  src_item = lst_string_get_item(src, src_index);
  dst_item = lst_string_get_item(dst, dst_index);

  char *csrc_item= (char*)src_item;
  char *cdst_item= (char*)dst_item;
  *cdst_item= *csrc_item;
//  src->sclass->copy_func(src_item, dst_item);
}


u_int
lst_string_items_common(LST_String *s1, u_int off1,
			LST_String *s2, u_int off2,
			u_int max_len)
{
  if (!s1 || !s2 || off1 >= s1->num_items || off2 >= s2->num_items)
	 return 0;

  const u_int len = MIN(s1->num_items-off1-1,
								MIN(s2->num_items-off2-1, max_len));

  u_int item1= off1;
  u_int item2= off2;
  for (u_int i = 0; i < len; ++i, ++item1, ++item2) {
	 if (!(*((char*)lst_string_get_item(s1, item1)) == *((char*)lst_string_get_item(s2, item2)))) {
		return i;
	 }
  }
  if (len==max_len)
	 return len;
  else {
	 if ((s1==s2) && (s1->num_items-1==item1) && (s2->num_items-1==item2))
		return len+1;
	 else
		return len;
  }
}


void
lst_string_index_init(LST_StringIndex *index)
{
  if (!index)
	 return;

  memset(index, 0, sizeof(LST_StringIndex));
  index->end_index = &index->end_index_local;
}


void
lst_string_index_copy(LST_StringIndex *src, LST_StringIndex *dst)
{
  if (!src || !dst)
	 return;

  dst->string = src->string;
  dst->start_index = src->start_index;
  *(dst->end_index) = *(src->end_index);
}


LST_StringSet *
lst_stringset_new(void)
{
  LST_StringSet *set;

  set = calloc(1, sizeof(LST_StringSet));

  if (!set)
	 return NULL;

  LIST_INIT(&set->members);

  return set;
}


void
lst_stringset_add(LST_StringSet *set, LST_String *string)
{
  if (!set || !string)
	 return;

  LIST_INSERT_HEAD(&set->members, string, set);
  set->size++;
}


void
lst_stringset_remove(LST_StringSet *set, LST_String *string)
{
  LST_String *set_string;

  if (!set || !string)
	 return;

  for (set_string = set->members.lh_first; set_string; set_string = set_string->set.le_next)
	 {
		if (set_string->id != string->id)
	continue;

		LIST_REMOVE(string, set);
		set->size--;
		return;
	 }
}


void           
lst_stringset_foreach(LST_StringSet *set,
		      LST_StringCB callback,
		      void *user_data)
{
  LST_String *string;

  if (!set || !callback)
    return;

  for (string = set->members.lh_first; string; string = string->set.le_next)
    callback(string, user_data);
}


void           
lst_stringset_free(LST_StringSet *set)
{
  LST_String *string;

  if (!set)
    return;

  while (set->members.lh_first)
    {
      string = set->members.lh_first;
      LIST_REMOVE(set->members.lh_first, set);
      lst_string_free(string);
    }

  free(set);
}
