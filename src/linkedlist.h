/* Header for linkedlist.c */

#ifndef LINKEDLIST_H_
#define LINKEDLIST_H_

// includes
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


 // structs
 typedef struct linkedlist_int_node linkedlist_int_node;
 struct linkedlist_int_node {
    int value;
    linkedlist_int_node *next;
};

// functions

// tests
void test_linkedlist(void);


#endif // LINKEDLIST_H_