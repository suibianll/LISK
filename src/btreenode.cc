
#include "btreenode.h"

namespace lisk{

    Node* BTree(_key_t *keys, _value_t *values,  uint64_t num, uint16_t depth){
        uint64_t num_leaves = (num + BTreeLeafNode::maxEntries - 1) / BTreeLeafNode::maxEntries;
        uint64_t num_items = num;

        std::pair<_key_t, uint64_t> *nextlevel = new std::pair<_key_t, uint64_t>[num_leaves];
        
        BTreeLeafNode* first_leaf = nullptr;
        BTreeLeafNode* last_leaf = nullptr;
        uint64_t cur = 0;
        for (size_t i = 0; i < num_leaves; ++i){
            BTreeLeafNode* leaf = new BTreeLeafNode(depth);
            uint64_t bulk_num = num_items /(num_leaves - i);
            for(size_t j = 0; j < bulk_num; ++j){
                leaf->append(keys[cur], values[cur]);
                ++cur;
            }
            nextlevel[i].first = keys[cur - 1]; 
            nextlevel[i].second = (uint64_t)leaf;
            if(last_leaf != nullptr){
                last_leaf->sibling_ = leaf;	
            }
            else
            first_leaf = leaf;
            last_leaf = leaf;

            num_items -= bulk_num; 
        }
        
        if(first_leaf == last_leaf){
            return (Node*) (SET_TYPE(first_leaf, BLEAF));
        }

        //build inner nodes
        uint64_t num_parents = num_leaves;
        while (num_parents != 1){
            uint64_t num_children = num_parents;
            num_parents = (num_children + BTreeInnerNode::maxEntries - 1) / (BTreeInnerNode::maxEntries - 1);
            uint64_t inner_index = 0;
            for (size_t i = 0; i < num_parents; i++){
            uint64_t bulk_num = num_children /(num_parents - i) - 1;
            BTreeInnerNode* new_inner = new BTreeInnerNode(depth);
            for (size_t j = 0; j < bulk_num; j++){
                new_inner->append(nextlevel[inner_index].first, (Node *)nextlevel[inner_index].second);
                ++inner_index;
            }
            new_inner->append((Node *)(nextlevel[inner_index].second));
            
            nextlevel[i].first = nextlevel[inner_index].first;
            nextlevel[i].second = (uint64_t)new_inner;
            ++inner_index;
            num_children -= bulk_num + 1;
            }
        }
        Node * root = (Node*) (SET_TYPE(nextlevel[0].second, BINNER));
        delete[] nextlevel;
        return root;
    }



    Node* makeRoot(_key_t k, Node* leftChild, Node* rightChild) {
          BTreeInnerNode* inner = new BTreeInnerNode(leftChild->depth);
          inner->key_num = 1;
          inner->datatslot[0].key = k;
          inner->datatslot[0].val = leftChild;
          inner->datatslot[1].val = rightChild;
          return (Node *) inner;
      }

    bool BTreeFind(Node* root, Key &key, _value_t &v, int depth){
        Node* node = root;
        _key_t k = key.getSlice(depth);
        // Parent of current node
        BTreeInnerNode* parent = nullptr;
        int search_depth = 0;
        while (node->node_type == BINNER) {
            BTreeInnerNode* inner = static_cast<BTreeInnerNode*>(node);
            node = inner->datatslot[inner->lowerBound(k)].val;
            node->Prefetch();
            #ifndef NDEBUG
            search_depth++;
            #endif
        }

        BTreeLeafNode* leaf = static_cast<BTreeLeafNode*>(node);
        leaf->Prefetch();
        unsigned pos = leaf->lowerBound(k);
        bool success = false;
        if ((pos < leaf->key_num) && (leaf->dataslot[pos].key==k)) {
            success = true;
            v = leaf->dataslot[pos].val;
        }
        #ifndef NDEBUG
        {
            index_depth += search_depth;
            btree_search ++;
        }
        #endif
        return success;
    }

    Node* BTreeInsert(Node* root, Key &key, _value_t v, int depth){
        // Current node
        Node* node = root;
        _key_t k = key.getSlice(depth);

        // Parent of current node
        BTreeInnerNode* parent = nullptr;
        Node* newRoot = nullptr;
        while (node->node_type==BINNER) {
            auto inner = static_cast<BTreeInnerNode *>(node);
            // Split eagerly if full
            if (inner->isFull()) {
                // Split
                _key_t sep; BTreeInnerNode* newInner = inner->split(sep);
                
                if (parent) parent->insert(sep , newInner);
                else    newRoot = makeRoot(sep, inner, newInner);
                if(k >= sep)	inner = newInner;
            }
            parent = inner;
            node = inner->datatslot[inner->lowerBound(k)].val;
          
        }

        BTreeLeafNode* leaf = static_cast<BTreeLeafNode*>(node);

        // Split leaf if full
        if (leaf->key_num == leaf->maxEntries) {
            // Split
            _key_t sep; BTreeLeafNode* newLeaf = leaf->split(sep);
            if(k > sep) newLeaf->insert(k, v);
            else leaf->insert(k, v);
                if (parent)
                    parent->insert(sep, newLeaf);
                else
                    newRoot = makeRoot(sep, leaf, newLeaf);
        } 
        else {
            leaf->insert(k, v);
            
        }
        return newRoot; // if split spreads to the root node, return the new root
    }

    bool BTreeUpdate(Node* root, Key &key, _value_t v, int depth){
        Node* node = root;
        _key_t k = key.getSlice(depth);
        // Parent of current node
        BTreeInnerNode* parent = nullptr;

        while (node->node_type == BINNER) {
            BTreeInnerNode* inner = static_cast<BTreeInnerNode*>(node);
            node = inner->datatslot[inner->lowerBound(k)].val;
        }

        BTreeLeafNode* leaf = static_cast<BTreeLeafNode*>(node);
        unsigned pos = leaf->lowerBound(k);
        bool success = false;
        if ((pos < leaf->key_num) && (leaf->dataslot[pos].key==k)) {
            success = true;
            leaf->dataslot[pos].val = v;
        }
        return success;
    }

    bool BTreeRemove(Node * root, Key &key, int depth){
        return false;
    }

    uint64_t BTreeScan(Key &key, int range, Record<_key_t, _value_t>* output, int depth){
        return 0;
    }


}