////////////////////////////////////////////////////////////////////////////
//
// * FILE: linkedlist.c
// * DESCRIPTION:
//    Linked lists.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 30/oct/2022
//
//////////////////////////////////////////////////////////////////////////

#include "linkedlist.h"



linkedlist_int_node *linkedlist_int_init(){
    linkedlist_int_node *head = malloc(sizeof(linkedlist_int_node));
    head->next = NULL;
    return head;
}

// Pushes new node with value `value` to the beginning of linked list pointed to by `*head`.
void linkedlist_int_push(linkedlist_int_node **head, int value){
    linkedlist_int_node *new_node = malloc(sizeof(linkedlist_int_node));
    new_node->value = value;
    new_node->next = *head;
    *head = new_node;
}

void test_linkedlist_int_push(void){
    linkedlist_int_node *head = NULL;
    linkedlist_int_push(&head, 2);
    linkedlist_int_push(&head, 10);
    assert((head->value == 10) && (head->next->value == 2));
    printf("Tests for `linkedlist_int_push` completed: no issues found.\n");
}

void test_linkedlist(void){
    test_linkedlist_int_push();
}