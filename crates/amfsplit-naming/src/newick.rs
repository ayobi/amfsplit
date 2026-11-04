//! Newick trees as Bio.Phylo 1.84 reads, prunes and writes them, for the withholding rounds: the frozen
//! external_holdout.py prunes the withheld tips from backbone.topology.newick with `Tree.prune` (in
//! sorted order) and writes the result with `format_branch_length='%1.15f'`.
use anyhow::{anyhow, bail, Result};

#[derive(Clone, Debug)]
struct Node { name: String, length: Option<f64>, children: Vec<usize> }

#[derive(Clone, Debug)]
pub struct Tree { nodes: Vec<Node>, root: usize }

impl Tree {
    pub fn parse(text: &str) -> Result<Self> {
        let b: Vec<u8> = text.bytes().filter(|c| !c.is_ascii_whitespace()).collect();
        let mut t = Tree { nodes: Vec::new(), root: 0 };
        let mut i = 0;
        t.root = t.clade(&b, &mut i)?;
        if b.get(i) != Some(&b';') { bail!("newick: expected ';' at byte {i}"); }
        Ok(t)
    }

    fn clade(&mut self, b: &[u8], i: &mut usize) -> Result<usize> {
        let id = self.nodes.len();
        self.nodes.push(Node { name: String::new(), length: None, children: Vec::new() });
        if b.get(*i) == Some(&b'(') {
            *i += 1;
            loop {
                let c = self.clade(b, i)?;
                self.nodes[id].children.push(c);
                match b.get(*i) {
                    Some(b',') => *i += 1,
                    Some(b')') => { *i += 1; break; }
                    _ => bail!("newick: expected ',' or ')' at byte {i}"),
                }
            }
        }
        let start = *i;
        while *i < b.len() && !matches!(b[*i], b':' | b',' | b')' | b';' | b'(') { *i += 1; }
        self.nodes[id].name = String::from_utf8(b[start..*i].to_vec())?;
        if b.get(*i) == Some(&b':') {
            *i += 1;
            let s = *i;
            while *i < b.len() && !matches!(b[*i], b',' | b')' | b';') { *i += 1; }
            self.nodes[id].length = Some(std::str::from_utf8(&b[s..*i])?.parse()?);
        }
        Ok(id)
    }

    pub fn tips(&self) -> Vec<String> {
        let mut out = Vec::new();
        let mut stack = vec![self.root];
        while let Some(n) = stack.pop() {
            if self.nodes[n].children.is_empty() { out.push(self.nodes[n].name.clone()); }
            for c in self.nodes[n].children.iter().rev() { stack.push(*c); }
        }
        out
    }

    /// Path from the root's child down to the named tip (Bio.Phylo `get_path`, root excluded).
    fn path(&self, tip: &str) -> Option<Vec<usize>> {
        fn walk(t: &Tree, n: usize, tip: &str, acc: &mut Vec<usize>) -> bool {
            for &c in &t.nodes[n].children {
                acc.push(c);
                if (t.nodes[c].children.is_empty() && t.nodes[c].name == tip) || walk(t, c, tip, acc) { return true; }
                acc.pop();
            }
            false
        }
        let mut acc = Vec::new();
        walk(self, self.root, tip, &mut acc).then_some(acc)
    }

    /// Bio.Phylo `TreeMixin.prune`: remove the tip; a parent left with one child is collapsed into it
    /// (the child's length gains the parent's), or at the root the child becomes the root, length None.
    pub fn prune(&mut self, tip: &str) -> Result<()> {
        let path = self.path(tip).ok_or_else(|| anyhow!("no tip {tip}"))?;
        let target = *path.last().unwrap();
        let parent = if path.len() == 1 { self.root } else { path[path.len() - 2] };
        let pos = self.nodes[parent].children.iter().position(|&c| c == target).unwrap();
        self.nodes[parent].children.remove(pos);
        if self.nodes[parent].children.len() == 1 {
            let child = self.nodes[parent].children[0];
            if parent == self.root {
                self.nodes[child].length = None;
                self.root = child;
            } else {
                if let Some(l) = self.nodes[child].length { self.nodes[child].length = Some(l + self.nodes[parent].length.unwrap_or(0.0)); }
                let grandparent = if path.len() < 3 { self.root } else { path[path.len() - 3] };
                let idx = self.nodes[grandparent].children.iter().position(|&c| c == parent).unwrap();
                self.nodes[grandparent].children[idx] = child;
            }
        }
        Ok(())
    }

    /// Bio.Phylo's Newick writer with `format_branch_length='%1.15f'`: every clade gets ':' and its length
    /// (a missing length written as 0), then ';' and a newline as `Phylo.write` ends the file.
    pub fn to_newick(&self) -> String {
        fn w(t: &Tree, n: usize, out: &mut String) {
            let node = &t.nodes[n];
            if !node.children.is_empty() {
                out.push('(');
                for (k, c) in node.children.iter().enumerate() { if k > 0 { out.push(','); } w(t, *c, out); }
                out.push(')');
            }
            out.push_str(&node.name);
            out.push_str(&format!(":{:.15}", node.length.unwrap_or(0.0)));
        }
        let mut s = String::new();
        w(self, self.root, &mut s);
        s.push_str(";\n");
        s
    }
}
