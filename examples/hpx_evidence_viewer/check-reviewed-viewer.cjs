const fs=require('node:fs'), vm=require('node:vm'), assert=require('node:assert/strict');
const html=fs.readFileSync(__dirname+'/reviewed-index.html','utf8');
const scripts=[...html.matchAll(/<script[^>]*>([\s\S]*?)<\/script>/gi)].map(x=>x[1]);
const elements=new Map(), events={};
function element(id){
 if(elements.has(id))return elements.get(id);
 const h=Number((html.match(new RegExp('id="'+id+'"[^>]*height="(\\d+)"'))||[])[1]||350);
 const ctx=new Proxy({}, {get:(o,k)=>o[k]??(()=>{}),set:(o,k,v)=>(o[k]=v,true)});
 const e={width:800,height:h,dataset:{},style:{},textContent:'',innerHTML:'',classList:{toggle(){}},listeners:{},
  getAttribute(k){return String(this[k]);},getBoundingClientRect(){return {width:800,height:parseFloat(this.style.height)||h,left:0,top:0};},
  getContext(){return ctx;}, addEventListener(k,f){(this.listeners[k]??=[]).push(f);}};
 elements.set(id,e);return e;
}
const sandbox={console,Math,JSON,window:{devicePixelRatio:2,addEventListener(k,f){(events[k]??=[]).push(f);}},document:{getElementById:element},requestAnimationFrame:f=>f(),fetch:()=>{throw Error('Network request blocked');}};
vm.createContext(sandbox);for(const script of scripts)vm.runInContext(script,sandbox);
for(const fn of events.load||[])fn();
assert.equal(element('inspWi').textContent,'NOT_AVAILABLE');
const canvas=element('cvCoord'), height=canvas.height;
for(let i=0;i<8;i++) for(const fn of events.resize||[])fn();
assert.equal(canvas.height,height);assert.equal(canvas.listeners.click.length,1);
canvas.listeners.click[0]({clientX:0});assert.equal(vm.runInContext('selectedResNum',sandbox),1);
canvas.listeners.click[0]({clientX:799.9});assert.equal(vm.runInContext('selectedResNum',sandbox),462);
vm.runInContext("setConstructMode('mature'); focus3DOnResidue(240)",sandbox);
assert.equal(vm.runInContext('getConstructPos(24)',sandbox),1);
const point=vm.runInContext('projectPoint(focusTarget.x,focusTarget.y,focusTarget.z,800,350)',sandbox);
assert.equal(point.px,400); assert.equal(point.py,175);
console.log('PASS: load, no-profile state, DPR resize stability, listener uniqueness, boundary selection, mature mapping and 3D focus (mock DOM; not full browser QA)');
