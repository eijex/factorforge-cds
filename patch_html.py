
path = 'C:/Work/eijex/factorforge/web/index.html'
with open(path, 'r', encoding='utf-8') as f:
    html = f.read()

target = '<div class=\
flex-grow
overflow-y-auto
p-6
space-y-8
custom-scrollbar
bg-slate-50/50
dark:bg-slate-900/50\>'

new_block = '''<div class=\flex-grow
overflow-y-auto
p-6
space-y-8
custom-scrollbar
bg-slate-50/50
dark:bg-slate-900/50\>
                <!-- Version 3.4.6 -->
                <div class=\relative
pl-8
border-l-2
border-emerald-500\>
                    <div
                        class=\absolute
-left-
[9px]
top-0
w-4
h-4
bg-emerald-500
rounded-full
border-4
border-white
dark:border-slate-900
shadow-[0_0_10px_rgba
16
185
129
0.5
]
\>
                    </div>
                    <div class=\flex
items-center
space-x-2
mb-2\>
                        <span
                            class=\px-2
py-0.5
bg-emerald-100
dark:bg-emerald-900/30
text-emerald-700
dark:text-emerald-400
text-[10px]
font-bold
rounded-md
uppercase\>Current</span>
                        <h3 class=\font-bold
text-slate-800
dark:text-white\>v3.4.6 &mdash; AgentOS &amp; SLLM Preview</h3>
                        <span class=\text-slate-400
text-[10px]\>2026-09-05</span>
                    </div>
                    <ul class=\text-xs
text-slate-600
dark:text-slate-400
space-y-2
list-disc
ml-4\>
                        <li><b>AgentOS Hard Gate:</b> The deterministic AgentOS evaluator is now officially exposed as a production-ready MCP tool.</li>
                        <li><b>FactorForge SLLM:</b> The 4.0 SLLM model is now available as a Research Preview.</li>
                    </ul>
                </div>
'''

html = html.replace(target, new_block)
with open(path, 'w', encoding='utf-8') as f:
    f.write(html)

