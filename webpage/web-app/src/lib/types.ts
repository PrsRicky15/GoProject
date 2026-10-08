export type Mode = 'line' | 'heatmap' | 'surface' | 'iso';
export type Axis = 'x' | 'y' | 'z';
export type FieldStyle = 'heatmap' | 'contour';
export type InputNumber = number | string;
export interface PlotDraft {
  mode: Mode;
  expression: string;
  resolution: InputNumber;
  bounds: Record<Axis, [InputNumber, InputNumber]>;
  level: InputNumber;
  style: FieldStyle;
}
export interface PlotRequest {
  mode: Mode;
  expression: string;
  resolution: number;
  bounds: Partial<Record<Axis, [number, number]>>;
  level: number;
  style: FieldStyle;
}
export interface PlotResult {
  settings: PlotRequest;
  coords: { x: Float64Array; y?: Float64Array; z?: Float64Array };
  values: Float64Array;
  count: number;
  invalid: number;
  min: number;
  max: number;
}
export interface RenderStatus { busy: boolean; error: string }
